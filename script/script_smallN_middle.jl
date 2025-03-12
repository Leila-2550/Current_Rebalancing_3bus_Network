using Pkg
using OpenDSSDirect
import LinearAlgebra: cis
using Plots
const _ODSS = OpenDSSDirect

# # Define OpenDSS circuit and model path
path = "/Users/raj055/Documents/GitHub/Current_rebalancing_onebus_full_network/data/Network_N_small" # network N is the original, copy 3 i chnged the load values

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
F_inv = [1 1 1; 1 α^2 α; 1 α α^2]

# fixing the tap position
tap_position = 1.0
OpenDSSDirect.Circuit.SetActiveElement("Transformer.tr1")
OpenDSSDirect.Transformers.Tap(tap_position)

OpenDSSDirect.Solution.Solve() # Solving the network
# voltages = _ODSS.Circuit.AllBusVMag()

######## Just for testing
all_elements = OpenDSSDirect.Circuit.AllElementNames()## returns all th elements names in the network
## To see each element is connected to what bus
element_name = "Transformer.tr1"  #  or Load.1 ## this way we see this element is connected to what bus and phase
_ODSS.Circuit.SetActiveElement(element_name)
connected_buses = _ODSS.CktElement.BusNames()
## Reading a load powers
_ODSS.Circuit.SetActiveElement("Load.1")
load68_powers = _ODSS.CktElement.Powers()
## ## Reading an element currents
_ODSS.Circuit.SetActiveElement("Transformer.tr1")
transformer_currents = _ODSS.CktElement.CurrentsMagAng()
#############

# Function to calculate the symmetrical components 
function calculate_symmetrical_components(phase_currents)
    return F * phase_currents
end

### 1- Measuring the voltage at the substation and compute the unbalance factor based on 2 approaches
_ODSS.Circuit.SetActiveBus("6687")#subsation bus
voltage_sub_before = _ODSS.Bus.PuVoltage()
 [round(v, digits=3) for v in voltage_sub_before] # rounding the measurement
voltage_sym_before = calculate_symmetrical_components(voltage_sub_before[1:3]) # symmetrical components of voltage
  [round(v, digits=4) for v in voltage_sym_before] # rounding the measurement

VUF1_before =(abs(voltage_sym_before[3])/abs(voltage_sym_before[2]))*100 # |v_negative|/|v_positive| 
VUF2_before=(sqrt(abs(voltage_sym_before[1]^2+ voltage_sym_before[3]^2)) / abs(voltage_sym_before[2]))*100 # sqrt(|v_zero|^2+|v_negative|^2)/|v_positive| 
##########
_ODSS.Circuit.SetActiveBus("sourcebus_11000")#subsation bus
voltage_sub_before = _ODSS.Bus.PuVoltage()

# 2- current at the substation
function get_substation_currents()
    _ODSS.Circuit.SetActiveElement("Transformer.tr1")
    currents = _ODSS.CktElement.Currents() # _ODSS.CktElement.CurrentsMagAng() # to get the magnitude and angle
    I_a = round(currents[5] , digits=2)
    I_b = round(currents[6] , digits=2)
    I_c = round(currents[7] , digits=2)
    I_n = round(currents[8] , digits=2)
    return [I_a, I_b, I_c, I_n]
end
I1_before= get_substation_currents()
I1_a_before, I1_b_before, I1_c_before, I1_n_before= I1_before
I1_sym_before = calculate_symmetrical_components(I1_before[1:3])# current I1 symmetrical components
[round(i, digits=3) for i in I1_sym_before]# rounding the measurement
## current unbalance factor
I_A = abs(I1_a_before)  # Magnitude of phase A
I_B = abs(I1_b_before)  # Magnitude of phase B
I_C = abs(I1_c_before)  # Magnitude of phase C
I1_avg = (I_A + I_B + I_C) / 3
PCBF_before = I1_avg / maximum([I_A, I_B, I_C]) #CUF= average current/max(Ia,Ib,Ic)

#### 3- calculating the load current I2
function calculate_I_line2()
    # Initialize total currents for each phase
    I2_a = 0.0 + 0.0im  # Phase A
    I2_b = 0.0 + 0.0im  # Phase B
    I2_c = 0.0 + 0.0im  # Phase C
    I2_n = 0.0 + 0.0im  #  neutral
    # List of lines connected to bus 6687
    lines = [
        "Line.line_75588639_6985_6732", #
    ]
    # Loop through each line
    for line in lines
        _ODSS.Circuit.SetActiveElement(line)
        currents = _ODSS.CktElement.Currents() #First four elements, currents flowing out from bus1 (6687) toward bus2 (6732). Last four elements represent the currents flowing into bus1 (6687) from bus2 (6732).
        #we considered the current from loads to substation
            I2_a = currents[1] # Phase A
            I2_b = currents[2] # Phase B
            I2_c = currents[3] # Phase C
            I2_n = currents[4] # Neutral
    end
    # Return the summed phase currents
    return [round(I2_a, digits=2), round(I2_b, digits=2), round(I2_c, digits=2), round(I2_n, digits=2)]
end
Iline2_before = calculate_I_line2()
Iline2_a_before, Iline2_b_before, Iline2_c_before,Iline2_n_before= Iline2_before
Iline2_0_before, Iline2_1_before,Iline2_2_before= calculate_symmetrical_components(Iline2_before[1:3])
##
function calculate_I2()
    # Initialize total currents for each phase
    I2_a = 0.0 + 0.0im  # Phase A
    I2_b = 0.0 + 0.0im  # Phase B
    I2_c = 0.0 + 0.0im  # Phase C
    I2_n = 0.0 + 0.0im  #  neutral
    # List of lines connected to bus 6687
    lines = [
        "Line.line_75588050_6732_6687", #
    ]
    # Loop through each line
    for line in lines
        _ODSS.Circuit.SetActiveElement(line)
        currents = _ODSS.CktElement.Currents() #First four elements, currents flowing out from bus1 (6687) toward bus2 (6732). Last four elements represent the currents flowing into bus1 (6687) from bus2 (6732).
        #we considered the current from loads to substation
            # 3-phase standard direction line
            I2_a += currents[1] # Phase A
            I2_b += currents[2] # Phase B
            I2_c += currents[3] # Phase C
            I2_n += currents[4] # Neutral
    
    end
    # Return the summed phase currents
    return [round(I2_a, digits=2), round(I2_b, digits=2), round(I2_c, digits=2), round(I2_n, digits=2)]
end
I2_before = calculate_I2()
I2_a_before, I2_b_before, I2_c_before,I2_n_before= I2_before
I2_0_before, I2_1_before,I2_2_before= calculate_symmetrical_components(I2_before[1:3])
#######3

#########
#get voltage and currents before compensation
#Retrieve all bus names in the network
all_buses = _ODSS.Circuit.AllBusNames()
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

# Get voltages before and after compensation
voltage_before = get_bus_voltages()
#################
# Define load elements
loads = ["Load.1", "Load.2", "Load.3", "Load.4", "Load.5", "Load.6"]

# Function to retrieve load currents (magnitude & phase)
function get_load_currents(loads)
    load_currents = Dict()
    load_phases = Dict()
    for load in loads
        _ODSS.Circuit.SetActiveElement(load)
        currents = _ODSS.CktElement.Currents()  # Complex currents
        current_mag_phase = [(abs(i), angle(i) * (180 / π)) for i in currents]  # Magnitude & phase (in degrees)
        # Store only the first current (assumed phase A)
        load_currents[load] = current_mag_phase[1][1]  # Magnitude
        load_phases[load] = current_mag_phase[1][2]    # Phase
    end
    return load_currents, load_phases
end
load_currents_before, load_phases_before = get_load_currents(loads)


## 4- Compensation Current I3
M1= [0 0 0; 0 0 0;0 0 1]#only zero sequence
M2= [1 0 0; 0 0 0;0 0 1]#negative and zero sequence ***

I3 = -F_inv * M2 * F * I2_before[1:3]
[round(i, digits=3) for i in I3]# rounding the measurement
I3_n=-(I3[1]+I3[2]+I3[3])
I3_sym = calculate_symmetrical_components(I3)
[round(i, digits=3) for i in I3_sym]# rounding the measurement
# Compute magnitude and phase angle (in degrees) for each phase
mag1_I3, mag2_I3, mag3_I3,= abs.(I3)  # Compute magnitude
ang1_I3, ang2_I3, ang3_I3= angle.(I3) .* (180 / π)  # Convert radians to degrees

### 5- adding compensator to the network as current source

# # Initialize empty arrays to store P and Q values for each phase
# S3 = zeros(ComplexF64, 3)  # Array to store complex apparent power
# P3 = zeros(Float64, 3)     # Array to store active power
# Q3 = zeros(Float64, 3)     # Array to store reactive power

OpenDSSDirect.Text.Command("New ISource.IDG1 Bus1=6732.1.4 Phases=1  , amps = $mag1_I3, ang = $ang1_I3")
OpenDSSDirect.Text.Command("New ISource.IDG2 Bus1=6732.2.4 Phases=1  , amps = $mag2_I3, ang = $ang2_I3")
OpenDSSDirect.Text.Command("New ISource.IDG3 Bus1=6732.3.4 Phases=1  , amps = $mag3_I3, ang = $ang3_I3")

all_elements = OpenDSSDirect.Circuit.AllElementNames()## returns all th elements names in the network

##7- solve the circuit and read the currents from current source to validate

OpenDSSDirect.Solution.Solve()
voltage_after = get_bus_voltages()
load_currents_after, load_phases_after = get_load_currents(loads)



## 5- Substation current after compensation
I1_after=[]
I1_after = -(I2_before[1:3]+I3)
[round(i, digits=3) for i in I1_after]# rounding the measurement
I1_a_after, I1_b_after , I1_c_after = I1_after
I1_n_after = -(I1_a_after+I1_b_after+I1_c_after)
I1_after_sym = calculate_symmetrical_components(I1_after)
ii = [round(i, digits=2) for i in I1_after_sym]

## current unbalance factor
I_A1 = abs(I1_a_after)  # Magnitude of phase A
I_B1 = abs(I1_b_after)  # Magnitude of phase B
I_C1 = abs(I1_c_after)  # Magnitude of phase C
I1_avg1 = (I_A1 + I_B1 + I_C1) / 3
CUF_after = I1_avg1 / maximum([I_A1, I_B1, I_C1]) #CUF= average current/max(Ia,Ib,Ic)


_ODSS.Circuit.SetActiveElement("ISource.IDG1")
current1 = _ODSS.CktElement.Currents()[2]
_ODSS.Circuit.SetActiveElement("ISource.IDG2")
current2 = OpenDSSDirect.CktElement.Currents()[2]
_ODSS.Circuit.SetActiveElement("ISource.IDG3")
current3 = OpenDSSDirect.CktElement.Currents()[2]

I3_measured = [current1, current2,current3]
[round(i, digits=3) for i in I3_measured]# rounding the measurement

_ODSS.Circuit.SetActiveBus("6687")
voltage_sub_after = _ODSS.Bus.PuVoltage() 
[round(v, digits=3) for v in voltage_sub_after]
# voltage_sub_after = measure_voltage("6687")
voltage_sym_aft = calculate_symmetrical_components(voltage_sub_after[1:3])
[round(v, digits=3) for v in voltage_sym_aft]
unbalance1_aft =abs(voltage_sym_aft[3])/abs(voltage_sym_aft[2])*100
unbalance2_aft=sqrt(abs(voltage_sym_aft[1]^2+ voltage_sym_aft[3]^2) / abs(voltage_sym_aft[2])^2 )*100

# Read the substation currents after compensation
_ODSS.Circuit.SetActiveElement("Transformer.tr1")
I1_trans_after1 = _ODSS.CktElement.Currents()[5:8]
I1_trans_after2 = get_substation_currents()

# Calculate symmetrical components of the compensated currents
I1_sym_trans_after = calculate_symmetrical_components(I1_trans_after2[1:3])
[round(i, digits=3) for i in I1_sym_trans_after]
# Check current unbalance factor (CUF) after compensation
I_A = abs(I1_trans_after2[1])
I_B = abs(I1_trans_after2[2])
I_C = abs(I1_trans_after2[3])
I1_avg_after = (I_A + I_B + I_C) / 3
CUF_after = I1_avg_after / maximum([I_A, I_B, I_C])
println("Current Unbalance Factor (CUF) After Compensation: ", CUF_after)



all_elements = OpenDSSDirect.Circuit.AllElementNames()## returns all th elements names in the network
print(all_elements)
_ODSS.Circuit.SetActiveElement("Generator.gen_a")
currents = OpenDSSDirect.CktElement.Currents()
println("Injected Currents: ", currents)

_ODSS.Circuit.SetActiveBus("6687")
voltage_sub_after = _ODSS.Bus.PuVoltage()
vv= [round(v, digits=3) for v in voltage_sub_after]
# voltage_sub_after = measure_voltage("6687")
voltage_sym_aft = calculate_symmetrical_components(voltage_sub_after[1:3])
unbalance2_aft=sqrt(abs(voltage_sym_aft[1]^2+ voltage_sym_aft[3]^2) / abs(voltage_sym_aft[2])^2 )*100
unbalance1_aft =abs(voltage_sym_aft[3])/abs(voltage_sym_aft[2])*100
########

###############
# Convert to arrays for plotting
bus_labels = collect(keys(voltage_before))
voltages_before = collect(values(voltage_before))
voltages_after = collect(values(voltage_after))


using Plots

# Convert Dict values to array for plotting
bus_labels = collect(keys(voltage_before))
voltages_before = [maximum(voltage_before[bus]) for bus in bus_labels]  # Take max magnitude for each bus
voltages_after = [maximum(voltage_after[bus]) for bus in bus_labels]

# --- 1. IMPROVED HISTOGRAM ---
histogram(voltages_before, bins=30, label="Before Compensation", alpha=0.5, color=:blue, xlabel="Voltage (p.u.)", ylabel="Frequency", title="Voltage Distribution Across Network")
histogram!(voltages_after, bins=30, label="After Compensation", alpha=0.5, color=:red)


## plotting voltages before and after compensation
# Ensure correct indexing for bus names
bus_indices = collect(1:length(bus_labels))  # Assign integer indices to buses

plot(bus_indices, voltages_before, 
    label="Before Compensation", 
    marker=:circle, 
    linestyle=:solid, 
    color=:blue, 
    xlabel="Bus", 
    ylabel="Voltage (p.u.)", 
    title="Voltage Profile Across Network", 
    legend=:topright, 
    xticks=(bus_indices, bus_labels),  # Ensure each bus gets correctly labeled
    rotation=45, 
    size=(600, 400)  # Wider figure for better readability
)

plot!(bus_indices, voltages_after, 
    label="After Compensation", 
    marker=:square, 
    linestyle=:dash, 
    color=:red
)
#####

# Convert Dict values to arrays for plotting
load_labels = collect(loads)  # Keep original load order
load_indices = collect(1:length(load_labels))  # Numeric indices to correctly place labels

currents_before_values = [load_currents_before[load] for load in load_labels]
currents_after_values = [load_currents_after[load] for load in load_labels]
phases_before_values = [load_phases_before[load] for load in load_labels]
phases_after_values = [load_phases_after[load] for load in load_labels]

# --- 1. CURRENT MAGNITUDE PROFILE ACROSS LOADS ---
plot(load_indices, currents_before_values, 
    label="Before Compensation", 
    marker=:circle, 
    linestyle=:solid, 
    color=:blue, 
    xlabel="Load", 
    ylabel="Current Magnitude (A)", 
    title="Load Current Profile", 
    legend=:topright, 
    xticks=(load_indices, load_labels),  # Ensure correct label assignment
    rotation=45, 
    size=(750,500)
)

plot!(load_indices, currents_after_values, 
    label="After Compensation", 
    marker=:square, 
    linestyle=:dash, 
    color=:red
)

# --- 2. CURRENT PHASE PROFILE ACROSS LOADS ---
plot(load_indices, phases_before_values, 
    label="Before Compensation", 
    marker=:circle, 
    linestyle=:solid, 
    color=:blue, 
    xlabel="Load", 
    ylabel="Phase Angle (Degrees)", 
    title="Load Current Phase Profile", 
    legend=:topright, 
    xticks=(load_indices, load_labels),  # Ensure correct label assignment
    rotation=45, 
    size=(750,500)
)

plot!(load_indices, phases_after_values, 
    label="After Compensation", 
    marker=:square, 
    linestyle=:dash, 
    color=:red
)
#####
###transformer current
# Extract magnitudes and phases
I1_after =[ -24.901 + 37.465im
44.896 + 2.832im
-19.995 - 40.297im
-0.00 + 0.00im]
I1_before_mags = [abs(i) for i in I1_before]
I1_after_mags = [abs(i) for i in I1_after]

I1_before_phases = [angle(i) * (180 / π) for i in I1_before]
I1_after_phases = [angle(i) * (180 / π) for i in I1_after]

# Define phase names
phases = ["Phase A", "Phase B", "Phase C", "Neutral"]

# --- 1. CURRENT MAGNITUDE COMPARISON ---
plot(phases, I1_before_mags, 
    label="Before Compensation", 
    marker=:circle, 
    linestyle=:solid, 
    color=:blue, 
    xlabel="Phase", 
    ylabel="Current Magnitude (A)", 
    title="Transformer Current Magnitude Before & After Compensation",
    legend=:topright
)

plot!(phases, I1_after_mags, 
    label="After Compensation", 
    marker=:square, 
    linestyle=:dash, 
    color=:red
)

### Plotting
using Unitful, Unitful.DefaultSymbols, PyPlot, ElectricalEngineering

a = 2.5  # Plot scale

# Constants for plotting
rc("text", usetex=false); rc("font", family="sans-serif", size=16)

# Plotting Substation Currents Before Compensation
println("Plotting Substation Currents Before Compensation")
phasorcosine(abs(I1_a_before), angle(I1_a_before), ylabel=L"$i$", maglabel=L"$\hat{I1}_a$", 
    labelrsep=0.5, figsize=(7*a, 2.5*a), color="blue", linestyle="-")
phasorcosine(abs(I1_b_before), angle(I1_b_before), ylabel=L"$i$", maglabel=L"$\hat{I1}_b$ ", 
    labelrsep=0.5, color="green", linestyle="-", add=true)
phasorcosine(abs(I1_c_before), angle(I1_c_before), ylabel=L"$i$", maglabel=L"$\hat{I1}_c$ ", 
    labelrsep=0.5, color="red", linestyle="-", add=true)
phasorcosine(abs(I1_n_before), angle(I1_n_before), ylabel=L"$i$", maglabel=L"$\hat{I1}_n$ ", 
    labelrsep=0.5, color="black", linestyle="--", add=true)
# Show the plot for Substation Currents Before Compensation
gcf()
show()

println("Plotting Substation Voltage before Compensation")
phasorcosine(abs(voltage_sub_before[1]), angle(voltage_sub_before[1]), ylabel=L"$i$", maglabel=L"$\hat{V1}_a$", 
    labelrsep=0.5, figsize=(7*a, 2.5*a), color="blue", linestyle="--")
phasorcosine(abs(voltage_sub_before[2]), angle(voltage_sub_before[2]), ylabel=L"$i$", maglabel=L"$\hat{V1}_b$", 
    labelrsep=0.5, color="green", linestyle="--", add=true)
phasorcosine(abs(voltage_sub_before[3]), angle(voltage_sub_before[3]), ylabel=L"$i$", maglabel=L"$\hat{V1}_c$ ", 
    labelrsep=0.5, color="red", linestyle="--", add=true)
phasorcosine(abs(voltage_sub_before[4]), angle(voltage_sub_before[4]), ylabel=L"$i$", maglabel=L"$\hat{V1}_n$ ", 
    labelrsep=0.5, color="red", linestyle="--", add=true)
# Show the plot for Substation Currents After Compensation
gcf()

# Plotting Substation Currents After Compensation
println("Plotting Substation Currents After Compensation")
phasorcosine(abs(I1_a_after), angle(I1_a_after), ylabel=L"$i$", maglabel=L"$\hat{I1}_a$", 
    labelrsep=0.5, figsize=(7*a, 2.5*a), color="blue", linestyle="--")
phasorcosine(abs(I1_b_after), angle(I1_b_after), ylabel=L"$i$", maglabel=L"$\hat{I1}_b$", 
    labelrsep=0.5, color="green", linestyle="--", add=true)
phasorcosine(abs(I1_c_after), angle(I1_c_after), ylabel=L"$i$", maglabel=L"$\hat{I1}_c$ ", 
    labelrsep=0.5, color="red", linestyle="--", add=true)
phasorcosine(abs(I1_n_after), angle(I1_n_after), ylabel=L"$i$", maglabel=L"$\hat{I1}_n$ ", 
    labelrsep=0.5, color="red", linestyle="--", add=true)
# Show the plot for Substation Currents After Compensation
gcf()
show()

println("Plotting Substation Voltage after Compensation")
phasorcosine(abs(voltage_sub_after[1]), angle(voltage_sub_after[1]), ylabel=L"$i$", maglabel=L"$\hat{V1}_a$", 
    labelrsep=0.5, figsize=(7*a, 2.5*a), color="blue", linestyle="--")
phasorcosine(abs(voltage_sub_after[2]), angle(voltage_sub_after[2]), ylabel=L"$i$", maglabel=L"$\hat{V1}_b$", 
    labelrsep=0.5, color="green", linestyle="--", add=true)
phasorcosine(abs(voltage_sub_after[3]), angle(voltage_sub_after[3]), ylabel=L"$i$", maglabel=L"$\hat{V1}_c$ ", 
    labelrsep=0.5, color="red", linestyle="--", add=true)
phasorcosine(abs(voltage_sub_after[4]), angle(voltage_sub_after[4]), ylabel=L"$i$", maglabel=L"$\hat{V1}_n$ ", 
    labelrsep=0.5, color="red", linestyle="--", add=true)
# Show the plot for Substation Currents After Compensation
gcf()

######

println("Plotting Substation Voltage before Compensation")
phasorcosine(abs(voltage_sub_before[1]), angle(voltage_sub_before[1]), ylabel=L"$i$", maglabel=L"$\hat{V1}_a$", 
    labelrsep=0.5, figsize=(7*a, 2.5*a), color="blue", linestyle="--")
phasorcosine(abs(voltage_sub_before[2]), angle(voltage_sub_before[2]), ylabel=L"$i$", maglabel=L"$\hat{V1}_b$", 
    labelrsep=0.5, color="green", linestyle="--", add=true)
phasorcosine(abs(voltage_sub_before[3]), angle(voltage_sub_before[3]), ylabel=L"$i$", maglabel=L"$\hat{V1}_c$ ", 
    labelrsep=0.5, color="red", linestyle="--", add=true)
phasorcosine(abs(voltage_sub_before[4]), angle(voltage_sub_before[4]), ylabel=L"$i$", maglabel=L"$\hat{V1}_n$ ", 
    labelrsep=0.5, color="red", linestyle="--", add=true)
# Show the plot for Substation Currents After Compensation
gcf()

println("Plotting Substation Voltage after Compensation")
phasorcosine(abs(voltage_sub_after[1]), angle(voltage_sub_after[1]), ylabel=L"$i$", maglabel=L"$\hat{V1}_a$", 
    labelrsep=0.5, figsize=(7*a, 2.5*a), color="blue", linestyle="--")
phasorcosine(abs(voltage_sub_after[2]), angle(voltage_sub_after[2]), ylabel=L"$i$", maglabel=L"$\hat{V1}_b$", 
    labelrsep=0.5, color="green", linestyle="--", add=true)
phasorcosine(abs(voltage_sub_after[3]), angle(voltage_sub_after[3]), ylabel=L"$i$", maglabel=L"$\hat{V1}_c$ ", 
    labelrsep=0.5, color="red", linestyle="--", add=true)
phasorcosine(abs(voltage_sub_after[4]), angle(voltage_sub_after[4]), ylabel=L"$i$", maglabel=L"$\hat{V1}_n$ ", 
    labelrsep=0.5, color="red", linestyle="--", add=true)
# Show the plot for Substation Currents After Compensation
gcf()



# Phasor diagram of substation currents before and after compensation
function plot_phasors(I_before, I_after)
    # plot(legend=:topright, title="Phasor Diagram of Substation Currents")
    plot!([0, real(I_before[1])], [0, imag(I_before[1])], arrow=:true, label="I_a Before", color=:blue)
    plot!([0, real(I_before[2])], [0, imag(I_before[2])], arrow=:true, label="I_b Before", color=:green)
    plot!([0, real(I_before[3])], [0, imag(I_before[3])], arrow=:true, label="I_c Before", color=:red)
    
    plot!([0, real(I_after[1])], [0, imag(I_after[1])], arrow=:true, label="I_a After", linestyle=:dash, color=:blue)
    plot!([0, real(I_after[2])], [0, imag(I_after[2])], arrow=:true, label="I_b After", linestyle=:dash, color=:green)
    plot!([0, real(I_after[3])], [0, imag(I_after[3])], arrow=:dash, label="I_c After", color=:red)
end
###########################################
