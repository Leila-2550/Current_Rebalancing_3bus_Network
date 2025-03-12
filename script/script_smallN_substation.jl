using Pkg
using OpenDSSDirect
import LinearAlgebra: cis
using Plots
const _ODSS = OpenDSSDirect

## Define OpenDSS circuit and model path
path = "/Users/raj055/Documents/GitHub/Current_Rebalancing_3bus_Network/data/Network_N_small" # network with 3buses:2 load bus and 1 substation bus

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

### 1- Measuring the voltage at the substation and compute the unbalance factor based on 2 approaches
_ODSS.Circuit.SetActiveBus("1")#subsation bus
U1_sub_before = _ODSS.Bus.PuVoltage()
 [round(v, digits=3) for v in U1_sub_before] # rounding the measurement
U1_sym_before = calculate_symmetrical_components(U1_sub_before[1:3]) # symmetrical components of voltage
  [round(v, digits=4) for v in U1_sym_before] # rounding the measurement

VUF1_before =(abs(U1_sym_before[3])/abs(U1_sym_before[2]))*100 # |v_negative|/|v_positive| 
VUF2_before=(sqrt(abs(U1_sym_before[1]^2+ U1_sym_before[3]^2)))/ abs(U1_sym_before[2])*100 # sqrt(|v_zero|^2+|v_negative|^2)/|v_positive| 

# voltage at bus 2
_ODSS.Circuit.SetActiveBus("2")#
U2_before = _ODSS.Bus.PuVoltage()
[round(v, digits=3) for v in U2_before]
U2_sym_before = calculate_symmetrical_components(U2_before[1:3]) # symmetrical components of voltage
  [round(v, digits=4) for v in U2_sym_before] # rounding the measurement
VUF1_2before =(abs(U2_sym_before[3])/abs(U2_sym_before[2]))*100 # |v_negative|/|v_positive| 
VUF2_2before=(sqrt(abs(U2_sym_before[1]^2+ U2_sym_before[3]^2)))/ abs(U2_sym_before[2])*100 # sqrt(|v_zero|^2+|v_negative|^2)/|v_positive| 

# voltage at bus 3
_ODSS.Circuit.SetActiveBus("3")
U3_before = _ODSS.Bus.PuVoltage()
[round(v, digits=3) for v in U3_before]
U3_sym_before = calculate_symmetrical_components(U3_before[1:3]) # symmetrical components of voltage
  [round(v, digits=4) for v in U2_sym_before] # rounding the measurement
VUF1_3before =(abs(U3_sym_before[3])/abs(U3_sym_before[2]))*100 # |v_negative|/|v_positive| 
VUF2_3before=(sqrt(abs(U3_sym_before[1]^2+ U3_sym_before[3]^2)))/ abs(U3_sym_before[2])*100 # sqrt(|v_zero|^2+|v_negative|^2)/|v_positive| 

##########
# 2- current at the substation and calculate the phase current balance factor (PCBF)
function get_substation_currents()
    _ODSS.Circuit.SetActiveElement("Transformer.tr1")
    currents = _ODSS.CktElement.Currents() # _ODSS.CktElement.CurrentsMagAng() # to get the magnitude and angle
    I_a = round(currents[5] , digits=2)
    I_b = round(currents[6] , digits=2)
    I_c = round(currents[7] , digits=2)
    I_n = round(currents[8] , digits=2)
    return [I_a, I_b, I_c, I_n]
end
I_sub_before= get_substation_currents()
I_sub_a_before, I_sub_b_before, I_sub_c_before,I_sub_n_before= I_sub_before
I_sub_sym_before = calculate_symmetrical_components(I_sub_before[1:3])# current I_sub symmetrical components
[round(i, digits=3) for i in I_sub_sym_before]# rounding the measurement
## current unbalance factor
I_A = abs(I_sub_a_before)  # Magnitude of phase A
I_B = abs(I_sub_b_before)  # Magnitude of phase B
I_C = abs(I_sub_c_before)  # Magnitude of phase C
I1_avg = (I_A + I_B + I_C) / 3
PCBF_before = I1_avg / maximum([I_A, I_B, I_C]) #CUF= average current/max(Ia,Ib,Ic)

#### 3- calculating the load current I_L12
function calculate_I2()
    _ODSS.Circuit.SetActiveElement("Line.line_75588050_2_1")
    currents = _ODSS.CktElement.Currents() # _ODSS.CktElement.CurrentsMagAng() # to get the magnitude and angle
    I_a = round(currents[5] , digits=2)
    I_b = round(currents[6] , digits=2)
    I_c = round(currents[7] , digits=2)
    I_n = round(currents[8] , digits=2)
    return [I_a, I_b, I_c, I_n]
end
I_l12_before = calculate_I2()
I_l12_0_before, I_l12_1_before,I_l12_2_before= calculate_symmetrical_components(I_l12_before[1:3])
#####
#current I_L23
function calculate_I3()
    _ODSS.Circuit.SetActiveElement("Line.line_75588639_3_2")
    currents = _ODSS.CktElement.Currents() # _ODSS.CktElement.CurrentsMagAng() # to get the magnitude and angle
    I_a = round(currents[5] , digits=2)
    I_b = round(currents[6] , digits=2)
    I_c = round(currents[7] , digits=2)
    I_n = round(currents[8] , digits=2)
    return [I_a, I_b, I_c, I_n]
end
I_l23_before = calculate_I3()
I_l23_0_before, I_l23_1_before,I_l23_2_before= calculate_symmetrical_components(I_l23_before[1:3])
#########
## 4- Calculating the Compensation Current I_comp
M1= [0 0 0; 0 0 0;0 0 1]# only zero sequence
M2= [1 0 0; 0 0 0;0 0 1]# negative and zero sequence ***

I_comp = -F_inv * M2 * F * I_l12_before[1:3]
[round(i, digits=3) for i in I_comp]# rounding the measurement
I_comp_n=-(I_comp[1]+I_comp[2]+I_comp[3])
I_comp_sym = calculate_symmetrical_components(I_comp)
[round(i, digits=3) for i in I_comp_sym]# rounding the measurement
# Compute magnitude and phase angle (in degrees) for each phase
mag1_I_comp, mag2_I_comp, mag3_I_comp= abs.(I_comp)  # Compute magnitude
ang1_I_comp, ang2_I_comp, ang3_I_comp= angle.(I_comp) .* (180 / π)  # Convert radians to degrees

########
## 5- Injecting the compensation current as current sources to the substation
OpenDSSDirect.Text.Command("New ISource.IDG1 Bus1=1.1.4 Phases=1  , amps = $mag1_I_comp, ang = $ang1_I_comp")
OpenDSSDirect.Text.Command("New ISource.IDG2 Bus1=1.2.4 Phases=1  , amps = $mag2_I_comp, ang = $ang2_I_comp")
OpenDSSDirect.Text.Command("New ISource.IDG3 Bus1=1.3.4 Phases=1  , amps = $mag3_I_comp, ang = $ang3_I_comp")

##validate that the current sources are added to the network
all_elements = OpenDSSDirect.Circuit.AllElementNames()## returns all th elements names in the network
########## 

#get voltage and currents before compensation before solving the network with the compensation
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

# Get bus voltage magnitudes before compensation
voltage_before = get_bus_voltages()

#############
## 6- solve the circuit and read the currents from current source to validate
OpenDSSDirect.Solution.Solve()

### 7- Substation current after compensation
I_sub_after = get_substation_currents()
# Calculate symmetrical components of the compensated currents
II_sub_sym_after = calculate_symmetrical_components(I_sub_after[1:3])
[round(i, digits=3) for i in II_sub_sym_after]
# Check PCBF after compensation
I_A = abs(I_sub_after[1])
I_B = abs(I_sub_after[2])
I_C = abs(I_sub_after[3])
I1_avg_after = (I_A + I_B + I_C) / 3
PCBF_after = I1_avg_after / maximum([I_A, I_B, I_C])
#currents at other lines
I_l12_after = calculate_I2()
I_l23_after = calculate_I3()

### 8- voltage at substation after
_ODSS.Circuit.SetActiveBus("1")
U1_sub_after = _ODSS.Bus.PuVoltage() 
[round(v, digits=3) for v in U1_sub_after]
# voltage_sub_after = measure_voltage("6687")
U1_sym_after = calculate_symmetrical_components(U1_sub_after[1:3])
[round(v, digits=3) for v in U1_sym_after]
VUF1_after = (abs(U1_sym_after[3])/abs(U1_sym_after[2]))*100
VUF2_after=(sqrt(abs(U1_sym_after[1]^2+ U1_sym_after[3]^2)))/ abs(U1_sym_after[2])*100 # sqrt(|v_zero|^2+|v_negative|^2)/|v_positive| 

# voltage at bus 2 after compensation
_ODSS.Circuit.SetActiveBus("2")#
U2_after = _ODSS.Bus.PuVoltage()
[round(v, digits=3) for v in U2_after]
U2_sym_after = calculate_symmetrical_components(U2_after[1:3])
[round(v, digits=3) for v in U2_sym_after]
VUF1_1after = (abs(U2_sym_after[3])/abs(U2_sym_after[2]))*100
VUF2_1after=(sqrt(abs(U2_sym_after[1]^2+ U2_sym_after[3]^2)))/ abs(U2_sym_after[2])*100 # sqrt(|v_zero|^2+|v_negative|^2)/|v_positive| 

# voltage at bus 3 after compensation
_ODSS.Circuit.SetActiveBus("3")
U3_after = _ODSS.Bus.PuVoltage()
[round(v, digits=3) for v in U3_after]
U3_sym_after = calculate_symmetrical_components(U3_after[1:3])
[round(v, digits=3) for v in U3_sym_after]
VUF1_2after = (abs(U3_sym_after[3])/abs(U3_sym_after[2]))*100
VUF2_2after=(sqrt(abs(U3_sym_after[1]^2+ U3_sym_after[3]^2)))/ abs(U3_sym_after[2])*100 # sqrt(|v_zero|^2+|v_negative|^2)/|v_positive| 

#### 
voltage_after = get_bus_voltages()#voltage magnitude at all buses after compensation
load_currents_after, load_phases_after = get_load_currents(loads) #currents at all loads after compensation


######## PLOTTING THE RESULTS #####
using Plots
##1-Plotting bus voltages per phase  before and after compensation
# Convert Dict values to arrays
bus_labels = collect(keys(voltage_before))  # Bus names
bus_indices = collect(1:length(bus_labels))  # Numeric indices for buses

# Extract phase voltages per bus
voltages_before_phase1 = [voltage_before[bus][1] for bus in bus_labels]  # Phase A
voltages_before_phase2 = [voltage_before[bus][2] for bus in bus_labels]  # Phase B
voltages_before_phase3 = [voltage_before[bus][3] for bus in bus_labels]  # Phase C

voltages_after_phase1 = [voltage_after[bus][1] for bus in bus_labels]  # Phase A
voltages_after_phase2 = [voltage_after[bus][2] for bus in bus_labels]  # Phase B
voltages_after_phase3 = [voltage_after[bus][3] for bus in bus_labels]  # Phase C

# LINE PLOT (Voltage Profile per Phase) ---
plot(bus_indices, voltages_before_phase1, label="Before - Phase A", marker=:circle, linestyle=:solid, color=:blue,
    xlabel="Bus", ylabel="Voltage (p.u.)", title="Voltage Profile (Before vs. After Compensation)", legend=:best,
    xticks=(bus_indices, bus_labels), rotation=45, size=(600, 400))

plot!(bus_indices, voltages_before_phase2, label="Before - Phase B", marker=:circle, linestyle=:solid, color=:green)
plot!(bus_indices, voltages_before_phase3, label="Before - Phase C", marker=:circle, linestyle=:solid, color=:red)

plot!(bus_indices, voltages_after_phase1, label="After - Phase A", marker=:square, linestyle=:dash, color=:lightblue)
plot!(bus_indices, voltages_after_phase2, label="After - Phase B", marker=:square, linestyle=:dash, color=:lightgreen)
plot!(bus_indices, voltages_after_phase3, label="After - Phase C", marker=:square, linestyle=:dash, color=:pink)

### 2- Plotting the line currents before and after compensation
# Define line names and phases
line_labels = ["L12-A", "L12-B", "L12-C", "L12-N",
               "L23-A", "L23-B", "L23-C", "L23-N"]

# Create a numeric index for each phase of each line
line_indices = collect(1:length(line_labels))


I_l12_before_mags = [abs(I) for I in I_l12_before]  # Magnitudes for line L12 before
I_l12_after_mags = [abs(I) for I in I_l12_after]    # Magnitudes for line L12 after

I_l23_before_mags = [abs(I) for I in I_l23_before]  # Magnitudes for line L23 before
I_l23_after_mags = [abs(I) for I in I_l23_after]    # Magnitudes for line L23 after
# Flatten current magnitudes into one list for proper plotting
I_before_mags = vcat(I_l12_before_mags, I_l23_before_mags)
I_after_mags  = vcat(I_l12_after_mags, I_l23_after_mags)


# --- 2.1. LINE PLOT (Line Current Profile per Phase) ---
plot(line_indices, I_before_mags, 
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
    size=(600, 500)
)

plot!(line_indices, I_after_mags, 
    label="After Compensation", 
    marker=:square, 
    linestyle=:dash, 
    color=:red
)


# Define phase labels
phase_labels = ["Phase A", "Phase B", "Phase C", "Neutral"]
phase_indices = collect(1:length(phase_labels))  # Numeric indices for phases

# --- 2.2 LINE PLOT for L12 ---
plot(phase_indices, I_l12_before_mags, 
    label="Before Compensation", 
    marker=:circle, 
    linestyle=:solid, 
    color=:blue,
    xlabel="Phase", 
    ylabel="Current Magnitude (A)", 
    title="Line L12 Current Profile (Before vs. After Compensation)", 
    legend=:best,
    xticks=(phase_indices, phase_labels), 
    size=(600, 400)
)

plot!(phase_indices, I_l12_after_mags, 
    label="After Compensation", 
    marker=:square, 
    linestyle=:dash, 
    color=:red
)

# --- 2.3 LINE PLOT for L23 ---
plot(phase_indices, I_l23_before_mags, 
    label="Before Compensation", 
    marker=:circle, 
    linestyle=:solid, 
    color=:blue,
    xlabel="Phase", 
    ylabel="Current Magnitude (A)", 
    title="Line L23 Current Profile (Before vs. After Compensation)", 
    legend=:best,
    xticks=(phase_indices, phase_labels), 
    size=(600, 400)
)

plot!(phase_indices, I_l23_after_mags, 
    label="After Compensation", 
    marker=:square, 
    linestyle=:dash, 
    color=:red
)

######
### 3- transformer current
I_sub_before_mags = [abs(i) for i in I_sub_before]
I_sub_after_mags = [abs(i) for i in I_sub_after]

I_sub_before_phases = [angle(i) * (180 / π) for i in I_sub_before]
I_sub_after_phases = [angle(i) * (180 / π) for i in I_sub_after]

# Define phase names
phases = ["Phase A", "Phase B", "Phase C", "Neutral"]

# --- 1. CURRENT MAGNITUDE COMPARISON ---
plot(phases, I_sub_before_mags, 
    label="Before Compensation", 
    marker=:circle, 
    linestyle=:solid, 
    color=:blue, 
    xlabel="Phase", 
    ylabel="Current Magnitude (A)", 
    title="Transformer Current Magnitude Before & After Compensation",
    legend=:best
)

plot!(phases, I_sub_after_mags, 
    label="After Compensation", 
    marker=:square, 
    linestyle=:dash, 
    color=:red
)

####4- Phasor plots
using Unitful, Unitful.DefaultSymbols, PyPlot, ElectricalEngineering
a = 2.5  # Plot scale
# Constants for plotting
rc("text", usetex=false); rc("font", family="sans-serif", size=16)

# Plotting Substation Currents Before Compensation
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

# Plotting Substation Currents after Compensation
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


println("Plotting Substation Voltage before Compensation")
phasorcosine(abs(U1_sub_before[1]), angle(U1_sub_before[1]), ylabel=L"$i$", maglabel=L"$\hat{V1}_a$", 
    labelrsep=0.5, figsize=(7*a, 2.5*a), color="blue", linestyle="--")
phasorcosine(abs(U1_sub_before[2]), angle(U1_sub_before[2]), ylabel=L"$i$", maglabel=L"$\hat{V1}_b$", 
    labelrsep=0.5, color="green", linestyle="--", add=true)
phasorcosine(abs(U1_sub_before[3]), angle(U1_sub_before[3]), ylabel=L"$i$", maglabel=L"$\hat{V1}_c$ ", 
    labelrsep=0.5, color="red", linestyle="--", add=true)
phasorcosine(abs(U1_sub_before[4]), angle(U1_sub_before[4]), ylabel=L"$i$", maglabel=L"$\hat{V1}_n$ ", 
    labelrsep=0.5, color="red", linestyle="--", add=true)

gcf()

println("Plotting Substation Voltage after Compensation")
phasorcosine(abs(U1_sub_after[1]), angle(U1_sub_after[1]), ylabel=L"$i$", maglabel=L"$\hat{V1}_a$", 
    labelrsep=0.5, figsize=(7*a, 2.5*a), color="blue", linestyle="--")
phasorcosine(abs(U1_sub_after[2]), angle(U1_sub_after[2]), ylabel=L"$i$", maglabel=L"$\hat{V1}_b$", 
    labelrsep=0.5, color="green", linestyle="--", add=true)
phasorcosine(abs(U1_sub_after[3]), angle(U1_sub_after[3]), ylabel=L"$i$", maglabel=L"$\hat{V1}_c$ ", 
    labelrsep=0.5, color="red", linestyle="--", add=true)
phasorcosine(abs(U1_sub_after[4]), angle(U1_sub_after[4]), ylabel=L"$i$", maglabel=L"$\hat{V1}_n$ ", 
    labelrsep=0.5, color="red", linestyle="--", add=true)

gcf()



println("Plotting  Voltage at bus 2 before Compensation")
phasorcosine(abs(U2_before[1]), angle(U2_before[1]), ylabel=L"$i$", maglabel=L"$\hat{V1}_a$", 
    labelrsep=0.5, figsize=(7*a, 2.5*a), color="blue", linestyle="--")
phasorcosine(abs(U2_before[2]), angle(U2_before[2]), ylabel=L"$i$", maglabel=L"$\hat{V1}_b$", 
    labelrsep=0.5, color="green", linestyle="--", add=true)
phasorcosine(abs(U2_before[3]), angle(U2_before[3]), ylabel=L"$i$", maglabel=L"$\hat{V1}_c$ ", 
    labelrsep=0.5, color="red", linestyle="--", add=true)
phasorcosine(abs(U2_before[4]), angle(U2_before[4]), ylabel=L"$i$", maglabel=L"$\hat{V1}_n$ ", 
    labelrsep=0.5, color="red", linestyle="--", add=true)

gcf()

println("Plotting Voltage at bus 2 after Compensation")
phasorcosine(abs(U2_after[1]), angle(U2_after[1]), ylabel=L"$i$", maglabel=L"$\hat{V1}_a$", 
    labelrsep=0.5, figsize=(7*a, 2.5*a), color="blue", linestyle="--")
phasorcosine(abs(U2_after[2]), angle(U2_after[2]), ylabel=L"$i$", maglabel=L"$\hat{V1}_b$", 
    labelrsep=0.5, color="green", linestyle="--", add=true)
phasorcosine(abs(U2_after[3]), angle(U2_after[3]), ylabel=L"$i$", maglabel=L"$\hat{V1}_c$ ", 
    labelrsep=0.5, color="red", linestyle="--", add=true)
phasorcosine(abs(U2_after[4]), angle(U2_after[4]), ylabel=L"$i$", maglabel=L"$\hat{V1}_n$ ", 
    labelrsep=0.5, color="red", linestyle="--", add=true)

gcf()