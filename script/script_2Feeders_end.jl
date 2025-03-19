using Pkg
using OpenDSSDirect
import LinearAlgebra: cis
using Plots
const _ODSS = OpenDSSDirect

## Define OpenDSS circuit and model path
path = "/Users/raj055/Documents/GitHub/Current_Rebalancing_3bus_Network/data/NetworkN_small_2Feeders" # network with 2 feeders

filename = path*"/Master.dss"
## Solving The Test Case (Unisng OpenDSS)
dss("""
    clear
    compile $filename
    solve
""")
dss("""
    clear
    compile $filename
    solve
""")

# Constants
α = cis(2π / 3)
F = (1 / 3) * [1 1 1; 1 α α^2; 1 α^2 α] # Fortescue transformation matrix
F_inv = [1 1 1; 1 α^2 α; 1 α α^2] # Inverse Fortescue transformation matrix

#  Tap position
tap_position = 1.0
OpenDSSDirect.Circuit.SetActiveElement("Transformer.tr1")
OpenDSSDirect.Transformers.Tap(tap_position)
OpenDSSDirect.Solution.Solve() # Solving the network

#####
# Function to calculate the symmetrical components 
function calculate_symmetrical_components(phase_currents)
    return F * phase_currents
end
# Function to calculate VUF
function calculate_VUF(U_sym)
    # Extract voltage sequence components
    v_zero = abs(U_sym[1])
    v_positive = abs(U_sym[2])
    v_negative = abs(U_sym[3])
    # Calculate Voltage Unbalance Factors
    VUF1 = (v_negative / v_positive) 
    VUF2 = (sqrt(v_zero^2 + v_negative^2) / v_positive) 
    return VUF1, VUF2
end
# Function to calculate PCBF
function calaculate_PCBF(I_value)
    I_A = abs(I_value[1])  # Mag of phase A
    I_B = abs(I_value[2])  # Mag of phase B
    I_C = abs(I_value[3])  # Mag of phase C
    I1_avg = (I_A + I_B + I_C) / 3
    PCBF= I1_avg / maximum([I_A, I_B, I_C]) #PCBF= average current/max(Ia,Ib,Ic)
    return PCBF
end

#function to measure the voltage at a specific bus
function measure_voltage(bus_name)
    _ODSS.Circuit.SetActiveBus(bus_name)  # Set the active bus
    voltage_magang = _ODSS.Bus.VMagAngle()  # Get voltage magnitudes and angles
    voltages= [round(v, digits=2) for v in voltage_magang] # rounding the measurement
    # Ensure the voltage vector has exactly 8 elements (3 phases + neutral)
    if length(voltage_magang) != 8
        error("Voltage vector must have exactly 8 elements (3 phases + neutral).")
    end
    va = (voltages[1],voltages[2]) # Phase A
    vb = (voltages[3],voltages[4]) # Phase B
    vc = (voltages[5],voltages[6]) # Phase C
    vn = (voltages[7],voltages[8]) # Neutral
    return va, vb, vc, vn
end

function convert_to_complex(input_tuple)
    value1 = [v[1] * cis(deg2rad(v[2])) for v in input_tuple]
    return [round(v, digits=2) for v in value1]
end

# Function to calculate the current at a specific element
function calculate_current(element_name)
    _ODSS.Circuit.SetActiveElement(element_name)  # Set the active element
    currents = _ODSS.CktElement.Currents()  # Get complex current values
    # Extract and round phase currents
    I_a = round(currents[5], digits=2)
    I_b = round(currents[6], digits=2)
    I_c = round(currents[7], digits=2)
    I_n = round(currents[8], digits=2)
    return [I_a, I_b, I_c, I_n] 
end

## function to calculate different voltgae parameters for a specific bus
function analyse_bus_voltage(bus_name)
    U = measure_voltage(bus_name)
    U_complex = convert_to_complex(U)
    U_symm = calculate_symmetrical_components(U_complex[1:3])
    U_sym = [round(v, digits=2) for v in U_symm] # rounding the measurement
    # voltage unbalance factor
    VUF = calculate_VUF(U_sym)
    VUF = [round(v, digits=4) for v in VUF] # rounding the measurement
    return U,U_complex, U_sym, VUF
end

## Function to analyse the current at a specific element
function analyse_current(element_name)
    # Get phase currents
    I = calculate_current(element_name)
    # Calculate symmetrical components
    I_symm = calculate_symmetrical_components(I[1:3])
    I_sym = [round(i, digits=2) for i in I_symm]  # Rounding
    # Calculate current unbalance factor
    PCBF = calaculate_PCBF(I)
    # Get current magnitudes and angles
    _ODSS.Circuit.SetActiveElement(element_name)
    I_mags = _ODSS.CktElement.CurrentsMagAng()[9:16]  # Extract magnitude and angle
    I_mags = [round(i, digits=2) for i in I_mags]  # Rounding
    return I,I_mags, I_sym, PCBF
end

# Function to retrieve voltage magnitudes for all buses
function get_bus_voltages()
    bus_voltages = Dict()
    for bus in all_buses
        _ODSS.Circuit.SetActiveBus(bus)
        voltage = _ODSS.Bus.PuVoltage()
        voltage_mag = [abs(v) for v in voltage[1:3]] # Get voltage magnitudes for three phases
        bus_voltages[bus] = voltage_mag  # Store all three phases instead of mean
    end
    return bus_voltages
end
##################
# Main Script
##################
##1- calculate voltage at all buses and calculate the voltage unbalance factor (VUF)
Bus_voltage_before = Dict()
#store voltage data in a dictionary
for bus in ["1", "2", "3", "4", "5"]
    U,U_complex, U_sym, VUF = analyse_bus_voltage(bus)
    Bus_voltage_before[bus] = (U=U,U_complex, Sym=U_sym, VUF=VUF)  # Store as named tuple
end
#printing results
for bus in ["1", "2", "3", "4", "5"]
    U,U_complex, U_sym, VUF = analyse_bus_voltage(bus)
    println("Bus $bus: U_before$bus:$U, Sym. Components_before$bus: $U_sym, VUF_before$bus: $VUF")
end

### 2- Calculate current and phase current balance factor (PCBF) at the substation and other lines 
Current_before = Dict()
#store current data in a dictionary
for element in ["Transformer.tr1", "Line.line_F1_2_1", "Line.line_F1_3_2", "Line.line_F2_4_1", "Line.line_F2_5_4"]
    I,I_mags, I_sym, PCBF= analyse_current(element)
    Current_before[element] = (I_complex= I, I_mag=I_mags, I_Sym=I_sym, PCBF=PCBF)  # Store as named tuple
end
#########

## 3- Calculating the Compensation Current I_comp
M1= [0 0 0; 0 0 0;0 0 1]# only zero sequence
M2= [1 0 0; 0 0 0;0 0 1]# negative and zero sequence ***

I_comp = -F_inv * M2 * F * (Current_before["Line.line_F1_2_1"].I_complex[1:3]+ Current_before["Line.line_F2_4_1"].I_complex[1:3])
[round(i, digits=2) for i in I_comp]# rounding the measurement
I_comp_n= sum(I_comp)
I_comp_sym = calculate_symmetrical_components(I_comp)
[round(i, digits=2) for i in I_comp_sym]# rounding the measurement
# Compute magnitude and phase angle (in degrees) for each phase
mag1_I_comp, mag2_I_comp, mag3_I_comp= abs.(I_comp)  # Compute magnitude
ang1_I_comp, ang2_I_comp, ang3_I_comp= angle.(I_comp) .* (180 / π)  # Convert radians to degrees

######## ADD COMPENSATOR ########
## 5- Injecting the compensation current as current sources to the end of feeder
# OpenDSSDirect.Text.Command("New ISource.IDG1 Bus1=3.1.4 Phases=1  , amps = $mag1_I_comp, ang = $ang1_I_comp")
# OpenDSSDirect.Text.Command("New ISource.IDG2 Bus1=3.2.4 Phases=1  , amps = $mag2_I_comp, ang = $ang2_I_comp")
# OpenDSSDirect.Text.Command("New ISource.IDG3 Bus1=3.3.4 Phases=1  , amps = $mag3_I_comp, ang = $ang3_I_comp")

mag_I_comp = abs.(I_comp)
ang_I_comp = angle.(I_comp) .* (180 / π)

### 4- Inject Compensation Current at the end of feeder 1
for (i, phase) in enumerate(["1.4", "2.4", "3.4"])
    OpenDSSDirect.Text.Command("New ISource.IDG$i Bus1=3.$phase Phases=1  , amps = $(mag_I_comp[i]), ang = $(ang_I_comp[i])")
end
##validate that the current sources are added to the network
all_elements = OpenDSSDirect.Circuit.AllElementNames()## returns all th elements names in the network
########## 

## 5- Get voltage magnitudes before compensation
all_buses = _ODSS.Circuit.AllBusNames()
# Get bus voltage magnitudes before compensation
voltage_before = get_bus_voltages()
###########

### 6- solve the circuit 
OpenDSSDirect.Solution.Solve()

### 7- Substation current after compensation
Current_after = Dict()
#store current data in a dictionary
for element in ["Transformer.tr1", "Line.line_F1_2_1", "Line.line_F1_3_2", "Line.line_F2_4_1", "Line.line_F2_5_4"]
    I,I_mags, I_sym, PCBF= analyse_current(element)
    Current_after[element] = (I_complex= I, I_mag=I_mags, I_Sym=I_sym, PCBF=PCBF)  # 
end
#printing results
for element in ["Transformer.tr1", "Line.line_F1_2_1", "Line.line_F1_3_2", "Line.line_F2_4_1", "Line.line_F2_5_4"]
    I,I_mags, I_sym, PCBF= analyse_current(element)
    println("Element$element: I_complex_after:$I, I_mag_after: $I_mags, I_sym_after: $I_sym, PCBF_after: $PCBF")
end
#### 8- voltage at all buses after compensation
Bus_voltage_after = Dict()
#store voltage data in a dictionary
for bus in ["1", "2", "3", "4", "5"]
    U,U_complex, U_sym, VUF = analyse_bus_voltage(bus)
    Bus_voltage_after[bus] = (U=U,U_complex, Sym=U_sym, VUF=VUF)  # 
end
#printing results
for bus in ["1", "2", "3", "4", "5"]
    U,U_complex, U_sym, VUF = analyse_bus_voltage(bus)
    println("Bus$bus: U_after:$U, Sym.Comp_after: $U_sym, VUF_after: $VUF")
end
## 9- Get voltage magnitudes after compensation
voltage_after = get_bus_voltages()


######## PLOTTING THE RESULTS #####
using Pkg
using Plots


## 1- Plotting Bus Voltages per Phase Before and After Compensation
bus_labels = ["sourcebus", "1", "2", "3", "4", "5"]
bus_indices = collect(1:length(bus_labels))  # Numeric indices for buses


# Extract phase voltages per bus
voltages_before_phase1 = [voltage_before[bus][1] for bus in bus_labels]  # Phase A
voltages_before_phase2 = [voltage_before[bus][2] for bus in bus_labels]  # Phase B
voltages_before_phase3 = [voltage_before[bus][3] for bus in bus_labels]  # Phase C

voltages_after_phase1 = [voltage_after[bus][1] for bus in bus_labels]  # Phase A
voltages_after_phase2 = [voltage_after[bus][2] for bus in bus_labels]  # Phase B
voltages_after_phase3 = [voltage_after[bus][3] for bus in bus_labels]  # Phase C

plot(bus_indices, voltages_before_phase1, label="Before - Phase A", marker=:circle, linestyle=:solid, color=:blue,
    xlabel="Bus", ylabel="Voltage (p.u.)", title="Voltage Magnitude Profile (Before vs. After Compensation)", legend=:bottomleft,
    xticks=(bus_indices, bus_labels), size=(600, 400))

plot!(bus_indices, voltages_before_phase2, label="Before - Phase B", marker=:circle, linestyle=:solid, color=:green)
plot!(bus_indices, voltages_before_phase3, label="Before - Phase C", marker=:circle, linestyle=:solid, color=:red)

plot!(bus_indices, voltages_after_phase1, label="After - Phase A", marker=:square, linestyle=:dash, color=:lightblue)
plot!(bus_indices, voltages_after_phase2, label="After - Phase B", marker=:square, linestyle=:dash, color=:lightgreen)
plot!(bus_indices, voltages_after_phase3, label="After - Phase C", marker=:square, linestyle=:dash, color=:pink)




using Plots

# Define the elements (lines and transformer) for which we have current data
elements = ["Transformer.tr1", "Line.line_F1_2_1", "Line.line_F1_3_2", "Line.line_F2_4_1", "Line.line_F2_5_4"]

# Extract current magnitudes before and after compensation
I_before = Dict(el => abs.(Current_before[el].I_complex) for el in elements)
I_after = Dict(el => abs.(Current_after[el].I_complex) for el in elements)

# --- 1. Plot Current Magnitude Comparison for Each Line ---

using Plots

# Ensure dictionaries contain data
if isempty(I_before) || isempty(I_after)
    error("Data in I_before or I_after is empty. Check the source dictionary.")
end

# Define labels for phases
phase_labels = ["Phase A", "Phase B", "Phase C", "Neutral"]
phase_indices = collect(1:length(phase_labels))

# Loop through each element and plot individual current magnitudes
for el in elements
    p = plot(phase_indices, I_before[el], 
        label="Before Compensation", 
        marker=:circle, 
        linestyle=:solid, 
        color=:blue,
        xlabel="Phase", 
        ylabel="Current Magnitude (A)", 
        title="Current Profile for $el (Before vs. After Compensation)", 
        legend=:best,
        xticks=(phase_indices, phase_labels), 
        size=(600, 400)
    )

    plot!(p, phase_indices, I_after[el], 
        label="After Compensation", 
        marker=:square, 
        linestyle=:dash, 
        color=:red
    )

    display(p)  # Show the plot
end


# --- 2. Combined Line Current Profile ---
line_labels = []
line_before_mags = []
line_after_mags = []

# Flatten data for combined plot
for el in elements[2:end]  # Exclude transformer
    for phase in 1:4
        push!(line_labels, "$el - $(phase_labels[phase])")
        push!(line_before_mags, I_before[el][phase])
        push!(line_after_mags, I_after[el][phase])
    end
end

line_indices = collect(1:length(line_labels))

plot(line_indices, line_before_mags, 
    label="Before Compensation", 
    marker=:circle, 
    linestyle=:solid, 
    color=:blue,
    xlabel="Line & Phase", 
    ylabel="Current Magnitude (A)", 
    title="Line Current Profile (Before vs. After Compensation)", 
    legend=:best,
    xticks=(line_indices, line_labels), 
    rotation=45, 
    size=(900, 500)
)

plot!(line_indices, line_after_mags, 
    label="After Compensation", 
    marker=:square, 
    linestyle=:dash, 
    color=:red
)

# --- 3. Transformer Current Magnitude Before and After Compensation ---
I_sub_before_mags = I_before["Transformer.tr1"]
I_sub_after_mags = I_after["Transformer.tr1"]

plot(phase_indices, I_sub_before_mags, 
    label="Before Compensation", 
    marker=:circle, 
    linestyle=:solid, 
    color=:blue, 
    xlabel="Phase", 
    ylabel="Current Magnitude (A)", 
    title="Transformer Current Magnitude Before & After Compensation",
    legend=:best
)

plot!(phase_indices, I_sub_after_mags, 
    label="After Compensation", 
    marker=:square, 
    linestyle=:dash, 
    color=:red
)

### ---4. Phasor plots
using Unitful, Unitful.DefaultSymbols, PyPlot, ElectricalEngineering
a = 2.5  # Plot scale
# Constants for plotting
rc("text", usetex=false); rc("font", family="sans-serif", size=16)

# Plotting Substation Currents Before Compensation
I_sub_before = Current_before["Transformer.tr1"].I_complex
println("Plotting Substation Currents Before Compensation")
phasorcosine(abs(I_sub_before[1]), angle(I_sub_before[1]), ylabel=L"$i$", maglabel=L"$\hat{I1}_a$", 
    labelrsep=0.5, figsize=(7*a, 2.5*a), color="blue", linestyle="-")
phasorcosine(abs(I_sub_before[2]), angle(I_sub_before[2]), ylabel=L"$i$", maglabel=L"$\hat{I1}_b$ ", 
    labelrsep=0.5, color="green", linestyle="-", add=true)
phasorcosine(abs(I_sub_before[3]), angle(I_sub_before[3]), ylabel=L"$i$", maglabel=L"$\hat{I1}_c$ ", 
    labelrsep=0.5, color="red", linestyle="-", add=true)
phasorcosine(abs(I_sub_before[4]), angle(I_sub_before[4]), ylabel=L"$i$", maglabel=L"$\hat{I1}_n$ ", 
    labelrsep=0.5, color="black", linestyle="--", add=true)
# Show the plot for Substation Currents Before Compensation
gcf()

# Plotting Substation Currents after Compensation
I_sub_after = Current_after["Transformer.tr1"].I_complex
println("Plotting Substation Currents after Compensation")
phasorcosine(abs(I_sub_after[1]), angle(I_sub_after[1]), ylabel=L"$i$", maglabel=L"$\hat{I1}_a$", 
    labelrsep=0.5, figsize=(7*a, 2.5*a), color="blue", linestyle="-")
phasorcosine(abs(I_sub_after[2]), angle(I_sub_after[2]), ylabel=L"$i$", maglabel=L"$\hat{I1}_b$ ", 
    labelrsep=0.5, color="green", linestyle="-", add=true)
phasorcosine(abs(I_sub_after[3]), angle(I_sub_after[3]), ylabel=L"$i$", maglabel=L"$\hat{I1}_c$ ", 
    labelrsep=0.5, color="red", linestyle="-", add=true)
phasorcosine(abs(I_sub_after[4]), angle(I_sub_after[4]), ylabel=L"$i$", maglabel=L"$\hat{I1}_n$ ", 
    labelrsep=0.5, color="black", linestyle="--", add=true)
gcf()


# Plotting  Currents at line 12 Before Compensation
I_sub_before = Current_before["Line.line_F1_2_1"].I_complex
println("Plotting Substation Currents Before Compensation")
phasorcosine(abs(I_sub_before[1]), angle(I_sub_before[1]), ylabel=L"$i$", maglabel=L"$\hat{I1}_a$", 
    labelrsep=0.5, figsize=(7*a, 2.5*a), color="blue", linestyle="-")
phasorcosine(abs(I_sub_before[2]), angle(I_sub_before[2]), ylabel=L"$i$", maglabel=L"$\hat{I1}_b$ ", 
    labelrsep=0.5, color="green", linestyle="-", add=true)
phasorcosine(abs(I_sub_before[3]), angle(I_sub_before[3]), ylabel=L"$i$", maglabel=L"$\hat{I1}_c$ ", 
    labelrsep=0.5, color="red", linestyle="-", add=true)
phasorcosine(abs(I_sub_before[4]), angle(I_sub_before[4]), ylabel=L"$i$", maglabel=L"$\hat{I1}_n$ ", 
    labelrsep=0.5, color="black", linestyle="--", add=true)
# Show the plot for Substation Currents Before Compensation
gcf()


# Plotting Currents  at line 12 after Compensation
I_sub_after = Current_after["Line.line_F1_2_1"].I_complex
println("Plotting Substation Currents after Compensation")
phasorcosine(abs(I_sub_after[1]), angle(I_sub_after[1]), ylabel=L"$i$", maglabel=L"$\hat{I1}_a$", 
    labelrsep=0.5, figsize=(7*a, 2.5*a), color="blue", linestyle="-")
phasorcosine(abs(I_sub_after[2]), angle(I_sub_after[2]), ylabel=L"$i$", maglabel=L"$\hat{I1}_b$ ", 
    labelrsep=0.5, color="green", linestyle="-", add=true)
phasorcosine(abs(I_sub_after[3]), angle(I_sub_after[3]), ylabel=L"$i$", maglabel=L"$\hat{I1}_c$ ", 
    labelrsep=0.5, color="red", linestyle="-", add=true)
phasorcosine(abs(I_sub_after[4]), angle(I_sub_after[4]), ylabel=L"$i$", maglabel=L"$\hat{I1}_n$ ", 
    labelrsep=0.5, color="black", linestyle="--", add=true)
gcf()


# Plotting  Currents at line 14 Before Compensation
I_sub_before = Current_before["Line.line_F2_4_1"].I_complex
println("Plotting Substation Currents Before Compensation")
phasorcosine(abs(I_sub_before[1]), angle(I_sub_before[1]), ylabel=L"$i$", maglabel=L"$\hat{I1}_a$", 
    labelrsep=0.5, figsize=(7*a, 2.5*a), color="blue", linestyle="-")
phasorcosine(abs(I_sub_before[2]), angle(I_sub_before[2]), ylabel=L"$i$", maglabel=L"$\hat{I1}_b$ ", 
    labelrsep=0.5, color="green", linestyle="-", add=true)
phasorcosine(abs(I_sub_before[3]), angle(I_sub_before[3]), ylabel=L"$i$", maglabel=L"$\hat{I1}_c$ ", 
    labelrsep=0.5, color="red", linestyle="-", add=true)
phasorcosine(abs(I_sub_before[4]), angle(I_sub_before[4]), ylabel=L"$i$", maglabel=L"$\hat{I1}_n$ ", 
    labelrsep=0.5, color="black", linestyle="--", add=true)
# Show the plot for Substation Currents Before Compensation
gcf()
show()

# Plotting Currents  at line 14 after Compensation
I_sub_after = Current_after["Line.line_F2_4_1"].I_complex
println("Plotting Substation Currents after Compensation")
phasorcosine(abs(I_sub_after[1]), angle(I_sub_after[1]), ylabel=L"$i$", maglabel=L"$\hat{I1}_a$", 
    labelrsep=0.5, figsize=(7*a, 2.5*a), color="blue", linestyle="-")
phasorcosine(abs(I_sub_after[2]), angle(I_sub_after[2]), ylabel=L"$i$", maglabel=L"$\hat{I1}_b$ ", 
    labelrsep=0.5, color="green", linestyle="-", add=true)
phasorcosine(abs(I_sub_after[3]), angle(I_sub_after[3]), ylabel=L"$i$", maglabel=L"$\hat{I1}_c$ ", 
    labelrsep=0.5, color="red", linestyle="-", add=true)
phasorcosine(abs(I_sub_after[4]), angle(I_sub_after[4]), ylabel=L"$i$", maglabel=L"$\hat{I1}_n$ ", 
    labelrsep=0.5, color="black", linestyle="--", add=true)
gcf()

#plotting the voltage at substation
U1_complex_before = Bus_voltage_before["1"].U_complex
println("Plotting Substation Voltage before Compensation")
phasorcosine(abs(U1_complex_before[1]), angle(U1_complex_before[1]), ylabel=L"$i$", maglabel=L"$\hat{U1}_a$", 
    labelrsep=0.5, figsize=(7*a, 2.5*a), color="blue", linestyle="-")
phasorcosine(abs(U1_complex_before[2]), angle(U1_complex_before[2]), ylabel=L"$i$", maglabel=L"$\hat{U1}_b$", 
    labelrsep=0.5, color="green", linestyle="-", add=true)
phasorcosine(abs(U1_complex_before[3]), angle(U1_complex_before[3]), ylabel=L"$i$", maglabel=L"$\hat{U1}_c$ ", 
    labelrsep=0.5, color="red", linestyle="-", add=true)
phasorcosine(abs(U1_complex_before[4]), angle(U1_complex_before[4]), ylabel=L"$i$", maglabel=L"$\hat{U1}_n$ ", 
    labelrsep=0.5, color="red", linestyle="-", add=true)
gcf()

U1_complex_aft = Bus_voltage_after["1"].U_complex
println("Plotting Substation Voltage after Compensation")
phasorcosine(abs(U1_complex_aft[1]), angle(U1_complex_aft[1]), ylabel=L"$i$", maglabel=L"$\hat{U1}_a$", 
    labelrsep=0.5, figsize=(7*a, 2.5*a), color="blue", linestyle="--")
phasorcosine(abs(U1_complex_aft[2]), angle(U1_complex_aft[2]), ylabel=L"$i$", maglabel=L"$\hat{U1}_b$", 
    labelrsep=0.5, color="green", linestyle="--", add=true)
phasorcosine(abs(U1_complex_aft[3]), angle(U1_complex_aft[3]), ylabel=L"$i$", maglabel=L"$\hat{U1}_c$ ", 
    labelrsep=0.5, color="red", linestyle="--", add=true)
phasorcosine(abs(U1_complex_aft[4]), angle(U1_complex_aft[4]), ylabel=L"$i$", maglabel=L"$\hat{U1}_n$ ", 
    labelrsep=0.5, color="red", linestyle="--", add=true)

gcf()

