using Pkg
using OpenDSSDirect
import LinearAlgebra: cis
using Plots
const _ODSS = OpenDSSDirect
using Statistics

# # Define OpenDSS circuit and model path
# path = "/Users/raj055/Documents/GitHub/Current_rebalancing_onebus_full_network/data/network_N" # network N is the original, copy 3 i chnged the load values
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
## Reading a bus voltage in p.u.
_ODSS.Circuit.SetActiveBus("6687")
voltage = _ODSS.Bus.PuVoltage()
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
U1_sub_before = _ODSS.Bus.PuVoltage()
[round(v, digits=3) for v in U1_sub_before] # rounding the measurement
U1_sym_before = calculate_symmetrical_components(U1_sub_before[1:3]) # symmetrical components of voltage
[round(v, digits=4) for v in U1_sym_before] # rounding the measurement

VUF1_before =(abs(U1_sym_before[3])/abs(U1_sym_before[2]))*100 # |v_negative|/|v_positive| 
VUF2_before=(sqrt(abs(U1_sym_before[1]^2+ U1_sym_before[3]^2)) / abs(U1_sym_before[2]))*100 # sqrt(|v_zero|^2+|v_negative|^2)/|v_positive| 
 
##########

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
I_sub_before= get_substation_currents()
I_sub_a_before, I_sub_b_before, I_sub_c_before, I_sub_n_before= I_sub_before
I_sub_sym_before = calculate_symmetrical_components(I_sub_before[1:3])# current I1 symmetrical components
[round(i, digits=3) for i in I_sub_sym_before]# rounding the measurement
## current unbalance factor
I_A = abs(I_sub_a_before)  # Magnitude of phase A
I_B = abs(I_sub_b_before)  # Magnitude of phase B
I_C = abs(I_sub_c_before)  # Magnitude of phase C
I1_avg = (I_A + I_B + I_C) / 3
PCBF_before = I1_avg / maximum([I_A, I_B, I_C]) #CUF= average current/max(Ia,Ib,Ic)

#### 3- calculating the load current I2
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
            I2_a = currents[5] # Phase A
            I2_b = currents[6] # Phase B
            I2_c = currents[7] # Phase C
            I2_n = currents[8] # Neutral
    
    end
    # Return the summed phase currents
    return [round(I2_a, digits=2), round(I2_b, digits=2), round(I2_c, digits=2), round(I2_n, digits=2)]
end
I_l12_before = calculate_I2()
I_l12_a_before, I_l12_b_before, I_l12_c_before,I_l12_n_before= I_l12_before
I_l12_0_before, I_l12_1_before,I_l12_2_before= calculate_symmetrical_components(I_l12_before[1:3])

_ODSS.Circuit.SetActiveElement("Line.line_75588639_6985_6732")
curr = _ODSS.CktElement.Currents()[5:8]
[round(i, digits=3) for i in curr]# rounding the measurement
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

I_comp = -F_inv * M2 * F * I_l12_before[1:3] #compensator current
[round(i, digits=3) for i in I_comp]# rounding the measurement
I_comp_n=-(I_comp[1]+I_comp[2]+I_comp[3])
I_comp_sym = calculate_symmetrical_components(I_comp)
[round(i, digits=3) for i in I_comp_sym]# rounding the measurement
# Compute magnitude and phase angle (in degrees) for each phase
mag1_I_comp, mag2_I_comp, mag3_I_comp,= abs.(I_comp)  # Compute magnitude
ang1_I_comp, ang2_I_comp, ang3_I_comp= angle.(I_comp) .* (180 / π)  # Convert radians to degrees



### 5- adding compensator to the network as current source

# # Initialize empty arrays to store P and Q values for each phase
# S3 = zeros(ComplexF64, 3)  # Array to store complex apparent power
# P3 = zeros(Float64, 3)     # Array to store active power
# Q3 = zeros(Float64, 3)     # Array to store reactive power

OpenDSSDirect.Text.Command("New ISource.IDG1 Bus1=6687.1.4 Phases=1  , amps = $mag1_I_comp, ang = $ang1_I_comp")
OpenDSSDirect.Text.Command("New ISource.IDG2 Bus1=6687.2.4 Phases=1  , amps = $mag2_I_comp, ang = $ang2_I_comp")
OpenDSSDirect.Text.Command("New ISource.IDG3 Bus1=6687.3.4 Phases=1  , amps = $mag3_I_comp, ang = $ang3_I_comp")

all_elements = OpenDSSDirect.Circuit.AllElementNames()## returns all th elements names in the network


##6- solve the circuit and read the currents from current source to validate

OpenDSSDirect.Solution.Solve()

#7- # Read the substation currents after compensation
_ODSS.Circuit.SetActiveElement("Transformer.tr1")
I1_trans_after1 = _ODSS.CktElement.Currents()[5:8]

I_sub_after = get_substation_currents()
# Calculate symmetrical components of the compensated currents
II_sub_sym_after = calculate_symmetrical_components(I_sub_after[1:3])
[round(i, digits=3) for i in II_sub_sym_after]
# Check current unbalance factor (CUF) after compensation
I_A = abs(I_sub_after[1])
I_B = abs(I_sub_after[2])
I_C = abs(I_sub_after[3])
I1_avg_after = (I_A + I_B + I_C) / 3
CUF_after = I1_avg_after / maximum([I_A, I_B, I_C])
println("Current Unbalance Factor (CUF) After Compensation: ", CUF_after)

I_l12_after = calculate_I2()


voltage_after = get_bus_voltages()
load_currents_after, load_phases_after = get_load_currents(loads)


_ODSS.Circuit.SetActiveElement("ISource.IDG1")
current1 = _ODSS.CktElement.Currents()[2]
_ODSS.Circuit.SetActiveElement("ISource.IDG2")
current2 = OpenDSSDirect.CktElement.Currents()[2]
_ODSS.Circuit.SetActiveElement("ISource.IDG3")
current3 = OpenDSSDirect.CktElement.Currents()[2]
I3_measured = [current1, current2,current3]


_ODSS.Circuit.SetActiveBus("6687")
U1_sub_after = _ODSS.Bus.PuVoltage() 
[round(v, digits=3) for v in U1_sub_after]
# voltage_sub_after = measure_voltage("6687")
U1_sub_sym_after = calculate_symmetrical_components(U1_sub_after[1:3])
[round(v, digits=3) for v in U1_sub_sym_after ]
unbalance1_aft =abs(U1_sub_sym_after[3])/abs(U1_sub_sym_after[2])*100
unbalance2_aft=sqrt(abs(U1_sub_sym_after[1]^2+ U1_sub_sym_after[3]^2) / abs(U1_sub_sym_after[2])^2 )*100

_ODSS.Circuit.SetActiveElement("Line.line_75588050_6732_6687")# Line.line_75588639_6985_6732
curr = _ODSS.CktElement.Currents()[5:8]
[round(i, digits=3) for i in curr]# rounding the measurement


all_elements = OpenDSSDirect.Circuit.AllElementNames()## returns all th elements names in the network
print(all_elements)
_ODSS.Circuit.SetActiveElement("Generator.gen_a")
currents = OpenDSSDirect.CktElement.Currents()
println("Injected Currents: ", currents)


_ODSS.Circuit.SetActiveBus("6732")
voltage_sub_after = _ODSS.Bus.PuVoltage()
vv= [round(v, digits=3) for v in voltage_sub_after]
# voltage_sub_after = measure_voltage("6687")
voltage_sym_aft = calculate_symmetrical_components(voltage_sub_after[1:3])
unbalance2_aft=sqrt(abs(voltage_sym_aft[1]^2+ voltage_sym_aft[3]^2) / abs(voltage_sym_aft[2])^2 )*100
unbalance1_aft =abs(voltage_sym_aft[3])/abs(voltage_sym_aft[2])*100


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
    legend=:best, 
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
    legend=:best, 
    xticks=(load_indices, load_labels),  # Ensure correct label assignment
    rotation=45, 
    size=(700,400)
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
I_sub_after =[-25.0 + 37.29im
44.95 + 2.74im
-19.99 - 40.3im
 0.04 + 0.26im
]
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


### Plotting
using Unitful, Unitful.DefaultSymbols, PyPlot, ElectricalEngineering

a = 2.5  # Plot scale

# Constants for plotting
rc("text", usetex=false); rc("font", family="sans-serif", size=16)

# Plotting Substation Currents Before Compensation
println("Plotting Substation Currents Before Compensation")
phasorcosine(abs(I_sub_a_before), angle(I_sub_a_before), ylabel=L"$i$", maglabel=L"$\hat{I1}_a$", 
    labelrsep=0.5, figsize=(7*a, 2.5*a), color="blue", linestyle="-")
phasorcosine(abs(I_sub_b_before), angle(I_sub_b_before), ylabel=L"$i$", maglabel=L"$\hat{I1}_b$ ", 
    labelrsep=0.5, color="green", linestyle="-", add=true)
phasorcosine(abs(I_sub_c_before), angle(I_sub_c_before), ylabel=L"$i$", maglabel=L"$\hat{I1}_c$ ", 
    labelrsep=0.5, color="red", linestyle="-", add=true)
phasorcosine(abs(I_sub_n_before), angle(I_sub_n_before), ylabel=L"$i$", maglabel=L"$\hat{I1}_n$ ", 
    labelrsep=0.5, color="black", linestyle="--", add=true)
# Show the plot for Substation Currents Before Compensation
gcf()
show()

# Plotting Substation Currents after Compensation
println("Plotting Substation Currents Before Compensation")
phasorcosine(abs(I_sub_after[1]), angle(I_sub_after[1]), ylabel=L"$i$", maglabel=L"$\hat{I1}_a$", 
    labelrsep=0.5, figsize=(7*a, 2.5*a), color="blue", linestyle="-")
phasorcosine(abs(I_sub_after[2]), angle(I_sub_after[2]), ylabel=L"$i$", maglabel=L"$\hat{I1}_b$ ", 
    labelrsep=0.5, color="green", linestyle="-", add=true)
phasorcosine(abs(I_sub_after[3]), angle(I_sub_after[3]), ylabel=L"$i$", maglabel=L"$\hat{I1}_c$ ", 
    labelrsep=0.5, color="red", linestyle="-", add=true)
phasorcosine(abs(I_sub_after[4]), angle(I_sub_after[4]), ylabel=L"$i$", maglabel=L"$\hat{I1}_n$ ", 
    labelrsep=0.5, color="black", linestyle="--", add=true)
# Show the plot for Substation Currents Before Compensation
gcf()
show()

println("Plotting Substation Voltage before Compensation")
phasorcosine(abs(U1_sub_before[1]), angle(U1_sub_before[1]), ylabel=L"$i$", maglabel=L"$\hat{V1}_a$", 
    labelrsep=0.5, figsize=(7*a, 2.5*a), color="blue", linestyle="--")
phasorcosine(abs(U1_sub_before[2]), angle(U1_sub_before[2]), ylabel=L"$i$", maglabel=L"$\hat{V1}_b$", 
    labelrsep=0.5, color="green", linestyle="--", add=true)
phasorcosine(abs(U1_sub_before[3]), angle(U1_sub_before[3]), ylabel=L"$i$", maglabel=L"$\hat{V1}_c$ ", 
    labelrsep=0.5, color="red", linestyle="--", add=true)
phasorcosine(abs(U1_sub_before[4]), angle(U1_sub_before[4]), ylabel=L"$i$", maglabel=L"$\hat{V1}_n$ ", 
    labelrsep=0.5, color="red", linestyle="--", add=true)
# Show the plot for Substation Currents After Compensation
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


# Plotting Load Currents (I2)
println("Plotting Load Currents (I2)")
phasorcosine(abs(I2_before[1]), angle(I2_before[1]), ylabel=L"$i$", maglabel=L"$\hat{I2}_a$", 
    labelrsep=0.5, figsize=(7*a, 2.5*a), color="blue", linestyle="-")
phasorcosine(abs(I2_before[2]), angle(I2_before[2]), ylabel=L"$i$", maglabel=L"$\hat{I2}_b$ ", 
    labelrsep=0.5, color="green", linestyle="-", add=true)
phasorcosine(abs(I2_before[3]), angle(I2_before[3]), ylabel=L"$i$", maglabel=L"$\hat{I2}_c$ ", 
    labelrsep=0.5, color="red", linestyle="-", add=true)
phasorcosine(abs(I2_n_before), angle(I2_n_before), ylabel=L"$i$", maglabel=L"$\hat{I2}_n$ ", 
    labelrsep=0.5, color="black", linestyle="--", add=true)
# Show the plot for Load Currents
gcf()
show()

println("Plotting Compensating Currents (I3)")
phasorcosine(abs(I3[1]), angle(I3[1]), ylabel=L"$i$", maglabel=L"$\hat{I3}_a$", 
    labelrsep=0.5, figsize=(7*a, 2.5*a), color="blue", linestyle="-")
phasorcosine(abs(I3[2]), angle(I3[2]), ylabel=L"$i$", maglabel=L"$\hat{I3}_b$", 
    labelrsep=0.5, color="green", linestyle="-", add=true)
phasorcosine(abs(I3[3]), angle(I3[3]), ylabel=L"$i$", maglabel=L"$\hat{I3}_c$ ", 
    labelrsep=0.5, color="red", linestyle="-", add=true)
phasorcosine(abs(I3_n), angle(I3_n), ylabel=L"$i$", maglabel=L"$\hat{I3}_c$ ", 
    labelrsep=0.5, color="black", linestyle="--", add=true)
# Show the plot for Compensating Currents
gcf()
show()








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
##compensating at the end of feeder


path = "/Users/raj055/Documents/GitHub/Current_rebalancing_onebus_full_network/data/Network_N" # network N is the original, copy 3 i chnged the load values
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
F = (1 / 3) * [1 1 1; 1 α α^2; 1 α^2 α]
F_inv = [1 1 1; 1 α^2 α; 1 α α^2]

tap_position = 1.0
OpenDSSDirect.Circuit.SetActiveElement("Transformer.tr1")
OpenDSSDirect.Transformers.Tap(tap_position)
OpenDSSDirect.Solution.Solve()
# voltages = _ODSS.Circuit.AllBusVMag()
println("Transformer tap set to: ", tap_position)


#############
##calculate the symmetrical components
function calculate_symmetrical_components(phase_currents)
    return F * phase_currents
end
#####voltage at the substation and unbalance factor
_ODSS.Circuit.SetActiveBus("6687")#subsation bus
voltage_sub_bef = _ODSS.Bus.PuVoltage()
v_sub= [round(v, digits=3) for v in voltage_sub_bef]
voltage_sym = calculate_symmetrical_components(voltage_sub_bef[1:3])
unbalance1_bef =abs(voltage_sym[3])/abs(voltage_sym[2])*100
unbalance2_bef=sqrt(abs(voltage_sym[1]^2+ voltage_sym[3]^2) / abs(voltage_sym[2])^2 )*100
##########
#current at the substation
function get_substation_currents()
    _ODSS.Circuit.SetActiveElement("Transformer.tr1")
    currents = _ODSS.CktElement.CurrentsMagAng()
    I_a = round(currents[1,5] * cis(deg2rad(currents[2,5])), digits=2)
    I_b = round(currents[1,6] * cis(deg2rad(currents[2,6])), digits=2)
    I_c = round(currents[1,7] * cis(deg2rad(currents[2,7])), digits=2)
    I_n = round(currents[1,8] * cis(deg2rad(currents[2,8])), digits=2)
    return [I_a, I_b, I_c, I_n]
end

I1_before= get_substation_currents()
I1_a_before, I1_b_before, I1_c_before, I1_n_before= I1_before


I1_sym_before = calculate_symmetrical_components(I1_before[1:3])

I_A = abs(I1_a_before)  # Magnitude of phase A
I_B = abs(I1_b_before)  # Magnitude of phase B
I_C = abs(I1_c_before)  # Magnitude of phase C
I1_avg = (I_A + I_B + I_C) / 3
PCBF = I1_avg / maximum([I_A, I_B, I_C])

I3=zeros(3,1)
AA = [I1_sym_before[1] 0 I1_sym_before[3]]
I3 = -F_inv * transpose(AA)

I3_n=-(I3[1]+I3[2]+I3[3])
I3_sym = calculate_symmetrical_components(I3)
############

# Compute I3 (real and imaginary components)
I3_a_real, I3_a_imag = real(I3[1]), imag(I3[1])
I3_b_real, I3_b_imag = real(I3[2]), imag(I3[2])
I3_c_real, I3_c_imag = real(I3[3]), imag(I3[3])
# Voltage magnitudes for each phase from v end of feeder
#####voltage at the substation and unbalance factor
_ODSS.Circuit.SetActiveBus("6694")#subsation bus
voltage_6694 = _ODSS.Bus.PuVoltage()
v_sub1= [round(v, digits=3) for v in voltage_6694]
voltage_sym = calculate_symmetrical_components(voltage_6694[1:3])
V_a = abs(v_sub1[1])
V_b = abs(v_sub1[2])
V_c = abs(v_sub1[3])
# Add compensator in OpenDSS
OpenDSSDirect.Text.Command("New Load.70 Bus1=6694.1.4 Phases=1  kV=0.2309 kW=$(I3_a_real * V_a) kvar=$(I3_a_imag * V_a)")
OpenDSSDirect.Text.Command("New Load.71 Bus1=6694.2.4 Phases=1  kV=0.2309 kW=$(I3_b_real * V_b) kvar=$(I3_b_imag * V_b)")
OpenDSSDirect.Text.Command("New Load.72 Bus1=6694.3.4 Phases=1  kV=0.2309 kW=$(I3_c_real * V_c) kvar=$(I3_c_imag * V_c)")

# Solve the circuit
OpenDSSDirect.Solution.Solve()

_ODSS.Circuit.SetActiveBus("6687")
voltage_sub_after = _ODSS.Bus.PuVoltage()
vv= [round(v, digits=3) for v in voltage_sub_after]
# voltage_sub_after = measure_voltage("6687")
voltage_sym_aft = calculate_symmetrical_components(voltage_sub_after[1:3])
unbalance2_aft=sqrt(abs(voltage_sym_aft[1]^2+ voltage_sym_aft[3]^2) / abs(voltage_sym_aft[2])^2 )*100
unbalance1_aft =abs(voltage_sym_aft[3])/abs(voltage_sym_aft[2])*100

_ODSS.Circuit.SetActiveElement("Transformer.tr1")
currents_sub_aft = _ODSS.CktElement.CurrentsMagAng()
I1_a_aft = round(currents_sub_aft[1,5] * cis(deg2rad(currents_sub_aft[2,5])), digits=2)
I1_b_aft = round(currents_sub_aft[1,6] * cis(deg2rad(currents_sub_aft[2,6])), digits=2)
I1_c_aft = round(currents_sub_aft[1,7] * cis(deg2rad(currents_sub_aft[2,7])), digits=2)
I1_n_aft = round(currents_sub_aft[1,8] * cis(deg2rad(currents_sub_aft[2,8])), digits=2)
I_sub_aft = [I1_a_aft, I1_b_aft, I1_c_aft, I1_n_aft]

I1_a_aft = I1_a_before+ I3[1]
I1_b_aft = I1_b_before+ I3[2]
I1_c_aft = I1_c_before+ I3[3]

I_sub_sym_aft = calculate_symmetrical_components(I_sub_aft[1:3])

I1_A = abs(I1_a_aft)  # Magnitude of phase A
I1_B = abs(I1_b_aft)  # Magnitude of phase B
I1_C = abs(I1_c_aft)  # Magnitude of phase C
I2_avg = (I1_A + I1_B + I1_C) / 3
PCBF = I2_avg / maximum([I1_A, I1_B, I1_C])


# Solve the circuit
OpenDSSDirect.Solution.Solve()

_ODSS.Circuit.SetActiveBus("6687")
voltage_sub_after = _ODSS.Bus.PuVoltage()
vv= [round(v, digits=3) for v in voltage_sub_after]
# voltage_sub_after = measure_voltage("6687")
voltage_sym_aft = calculate_symmetrical_components(voltage_sub_after[1:3])
unbalance2_aft=sqrt(abs(voltage_sym_aft[1]^2+ voltage_sym_aft[3]^2) / abs(voltage_sym_aft[2])^2 )*100
unbalance1_aft =abs(voltage_sym_aft[3])/abs(voltage_sym_aft[2])*100








################################################



##compensating loads
function measure_voltage(bus_name)
    _ODSS.Circuit.SetActiveBus(bus_name)
    voltages_pu = _ODSS.Bus.PuVoltage()
    voltages = [(abs(v), angle(v) * 180 / π) for v in voltages_pu]
    return voltages
end
function apply_compensation(bus_name, I3, voltage_at_bus)
    # Initialize compensating powers for each phase
    compensating_loads = Dict()
    # Calculate P and Q for each phase
    for phase in 1:3
        # Extract voltage magnitude and angle
        V_mag = voltage_at_bus[phase][1]  # Magnitude of voltage
        V_ang = deg2rad(voltage_at_bus[phase][2])  # Angle in radians
        # Convert voltage to complex form
        V = V_mag * cis(V_ang)
        # Compensating current for this phase
        I = I3[phase]
        # Calculate complex power
        S = V * conj(I)  # Complex power = V * I_conjugate
        P = real(S)  # Active power
        Q = imag(S)  # Reactive power
        # Store compensating loads for each phase
        compensating_loads[phase] = (P, Q)
    end
    return compensating_loads
end

#####

abs(I1_sym_before[1])
angle(I1_sym_before[1]) * 180 / π
#### Method 2: calculate the symmetrical components
function calculate_symmetrical_components(phase_currents)
    return F * phase_currents
end
I1_sym_before = calculate_symmetrical_components(I1_before[1:3])
I_comp2 = I1_sym_before[1] + I1_sym_before[3]

abs(I1_sym_before[1])
angle(I1_sym_before[1]) * 180 / π

#####voltage at the substation and unbalance factor
_ODSS.Circuit.SetActiveBus("6687")#subsation bus
voltage_sub_bef = _ODSS.Bus.PuVoltage()
vr= [round(v, digits=3) for v in voltage_sub_bef]
voltage_sym = calculate_symmetrical_components(voltage_sub_bef[1:3])
unbalance1_bef =abs(voltage_sym[3])/abs(voltage_sym[2])*100
unbalance2_bef=sqrt(abs(voltage_sym[1]^2+ voltage_sym[3]^2) / abs(voltage_sym[2])^2 )*100
##########
function get_bus_voltage_snap()
    bus_dict = Dict()
    phase_mapping = Dict(1 => "a", 2 => "b", 3 => "c", 4 => "n")

    for bus_name in OpenDSSDirect.Circuit.AllBusNames()
        bus_dict[bus_name] = Dict()
        OpenDSSDirect.Circuit.SetActiveBus(bus_name)

        bus_phases = OpenDSSDirect.Bus.Nodes()  # Get the connected phases of the bus
        bus_voltages = OpenDSSDirect.Bus.PuVoltage()  # Get per-unit voltages of the bus
        # Check if neutral (phase 4) is present
        has_neutral = 4 ∈ bus_phases
        if has_neutral
            neutral_voltage = bus_voltages[length(bus_voltages)] 
        else
            neutral_voltage = 0+0im
        end
        # Remove phase 4 from the list of bus phases for processing
        bus_phases = filter(x -> x != 4, bus_phases)
        for (i, phase) in enumerate(bus_phases)
            voltage_phase = bus_voltages[i]
            # Calculate voltage magnitude and angle, with or without neutral reference
            if has_neutral
                bus_dict[bus_name]["vm" * phase_mapping[phase]] = abs(voltage_phase - neutral_voltage)
                bus_dict[bus_name]["va" * phase_mapping[phase]] = angle(voltage_phase - neutral_voltage) * 180 / pi
            else
                bus_dict[bus_name]["vm" * phase_mapping[phase]] = abs(voltage_phase)
                bus_dict[bus_name]["va" * phase_mapping[phase]] = angle(voltage_phase) * 180 / pi
            end
        end
    end
    return bus_dict
end
# Measure system voltages
bus_dict1 = get_bus_voltage_snap()  # Get the voltages of all phases of all buses
####################
function calculate_I2()
    # Initialize total currents for each phase
    I2_a = 0.0 + 0.0im  # Phase A
    I2_b = 0.0 + 0.0im  # Phase B
    I2_c = 0.0 + 0.0im  # Phase C
    # List of lines connected to bus 6687
    lines = [
        "Line.line_75588050_6732_6687",
        "Line.line_75588636_7025_6687",
        "Line.line_75588890_7691_6687",
        "Line.line_75589062_7785_6687",
        "Line.line_75589080_6687_7153"  # Reversed direction
    ]
    # Loop through each line
    for line in lines
        _ODSS.Circuit.SetActiveElement(line)
        currents = _ODSS.CktElement.CurrentsMagAng()
        if line == "Line.line_75589080_6687_7153"
            # Reversed direction line: use second half of the current matrix
            I2_a += currents[1, 5] * cis(deg2rad(currents[2, 5]))  # Phase A
            I2_b += currents[1, 6] * cis(deg2rad(currents[2, 6]))  # Phase B
            I2_c += currents[1, 7] * cis(deg2rad(currents[2, 7]))  # Phase C
        elseif line == "Line.line_75588636_7025_6687"
            # Phase A only (standard direction)
            I2_a -= currents[1, 1] * cis(deg2rad(currents[2, 1]))
        elseif line == "Line.line_75588890_7691_6687"###########load connected
            # Phase B only (standard direction)
            I2_b -= currents[1, 1] * cis(deg2rad(currents[2, 1]))
        elseif line == "Line.line_75589062_7785_6687"########### load connected
            # Phase C only (standard direction)
            I2_c -= currents[1, 1] * cis(deg2rad(currents[2, 1]))
        else
            # 3-phase standard direction line
            I2_a -= currents[1, 1] * cis(deg2rad(currents[2, 1]))  # Phase A
            I2_b -= currents[1, 2] * cis(deg2rad(currents[2, 2]))  # Phase B
            I2_c -= currents[1, 3] * cis(deg2rad(currents[2, 3]))  # Phase C
        end
    end

    # Return the summed phase currents
    return [round(I2_a, digits=2), round(I2_b, digits=2), round(I2_c, digits=2)]
end

I2_before = calculate_I2()
I2_n_before =round(-(I2_before[1]+I2_before[2]+I2_before[3]),digits=2) 

I2_0_before, I2_1_before,I2_2_before= calculate_symmetrical_components(I2_before)

######### Compensating current

M1= [0 0 0; 0 0 0;0 0 1]#only zero sequence
M2= [1 0 0; 0 0 0;0 0 1]#negative and zero sequence

I3 = -F_inv * M1 * F * I2_before
I3_n=-(I3[1]+I3[2]+I3[3])
I3_sym = calculate_symmetrical_components(I3)

I1_after=[]
I1_after = -(I2_before+I3)
I1_a_after, I1_b_after , I1_c_after = I1_after
I1_n_after = -(I1_a_after+I1_b_after+I1_c_after)
I1_after_sym = calculate_symmetrical_components(I1_after)
ii = [round(i, digits=2) for i in I1_after_sym]

### copensating load
function measure_voltage(bus_name)
    _ODSS.Circuit.SetActiveBus(bus_name)
    voltages_pu = _ODSS.Bus.PuVoltage()
    voltages = [(abs(v), angle(v) * 180 / π) for v in voltages_pu]
    return voltages
end
function apply_compensation(bus_name, I3, voltage_at_bus)
    # Initialize compensating powers for each phase
    compensating_loads = Dict()
    # Calculate P and Q for each phase
    for phase in 1:3
        # Extract voltage magnitude and angle
        V_mag = voltage_at_bus[phase][1]  # Magnitude of voltage
        V_ang = deg2rad(voltage_at_bus[phase][2])  # Angle in radians
        # Convert voltage to complex form
        V = V_mag * cis(V_ang)
        # Compensating current for this phase
        I = I3[phase]
        # Calculate complex power
        S = V * conj(I)  # Complex power = V * I_conjugate
        P = real(S)  # Active power
        Q = imag(S)  # Reactive power
        # Store compensating loads for each phase
        compensating_loads[phase] = (P, Q)
    end
    return compensating_loads
end
# Measure voltage at bus 6687
voltage_at_6687 = measure_voltage("6687")
# Apply compensation using I3
compensating_loads = apply_compensation("6687", I3, voltage_at_6687)

# Update compensating loads in OpenDSS
_ODSS.Circuit.SetActiveElement("70")  # Phase 1
_ODSS.Loads.Name("70")  # Phase 1
_ODSS.Loads.kW(compensating_loads[3][1])
_ODSS.Loads.kvar(compensating_loads[3][2])

_ODSS.Circuit.SetActiveElement("71")  # Phase 2
_ODSS.Loads.Name("71")  # Phase 2
_ODSS.Loads.kW(compensating_loads[1][1])
_ODSS.Loads.kvar(compensating_loads[1][2])

_ODSS.Circuit.SetActiveElement("72")  # Phase 3
_ODSS.Loads.Name("72")  # Phase 3
_ODSS.Loads.kW(compensating_loads[2][1])
_ODSS.Loads.kvar(compensating_loads[2][2])

# Solve the circuit again after applying compensation
_ODSS.Solution.Solve()

bus_dict2 = get_bus_voltage_snap()
_ODSS.Circuit.SetActiveBus("6687")
voltage_sub_after = _ODSS.Bus.PuVoltage()
vv= [round(v, digits=3) for v in voltage_sub_after]
# voltage_sub_after = measure_voltage("6687")
voltage_sym_aft = calculate_symmetrical_components(voltage_sub_after[1:3])
unbalance2_aft=sqrt(abs(voltage_sym_aft[1]^2+ voltage_sym_aft[3]^2) / abs(voltage_sym_aft[2])^2 )*100
unbalance1_aft =abs(voltage_sym_aft[3])/abs(voltage_sym_aft[2])*100

_ODSS.Circuit.SetActiveElement("Transformer.tr1")
currents_sub_aft = _ODSS.CktElement.CurrentsMagAng()
I_a_aft = round(currents_sub_aft[1,5] * cis(deg2rad(currents_sub_aft[2,5])), digits=2)
I_b_aft = round(currents_sub_aft[1,6] * cis(deg2rad(currents_sub_aft[2,6])), digits=2)
I_c_aft = round(currents_sub_aft[1,7] * cis(deg2rad(currents_sub_aft[2,7])), digits=2)
I_n_aft = round(currents_sub_aft[1,8] * cis(deg2rad(currents_sub_aft[2,8])), digits=2)
I_sub_aft = [I_a_aft, I_b_aft, I_c_aft, I_n_aft]

I_sub_sym_aft = calculate_symmetrical_components(I_sub_aft[1:3])
unbalance1_current_aft =abs(I_sub_sym_aft[3])/abs(I_sub_sym_aft[2])*100
unbalance2_current_aft=sqrt(abs(I_sub_sym_aft[1]^2+ I_sub_sym_aft[3]^2) / abs(I_sub_sym_aft[2])^2 )*100

# function inject_compensating_current(I3) ## calculate the compensating loads
#     # Define the names of the adjustable compensator loads
#     compensator_loads = ["CompLoad1", "CompLoad2", "CompLoad3"]
#     phases = ["1.1", "1.2", "1.3"]  # Phases a, b, and c at the slack bus

#     # Loop through each phase and adjust the corresponding compensator load
#     for i in 1:3
#         # Get the compensating current for the current phase (I3 is complex)
#         I3_phase = I3[i]
#         # Measure the voltage at the slack bus for the current phase
#         _ODSS.Circuit.SetActiveBus("6687")
#         voltage = _ODSS.Bus.PuVoltage()[1] * 0.23  # Convert pu voltage to kV (0.4 kV base)

#         # Calculate the required active (kW) and reactive (kvar) power
#         P_comp = voltage * real(I3_phase)
#         Q_comp = voltage * imag(I3_phase)

#         # Extract only the real part and ensure they are Float64
#         P_comp = Float64(real(P_comp))
#         Q_comp = Float64(real(Q_comp))

#         # Set the active and reactive power for the compensator load
#         _ODSS.Circuit.SetActiveElement("Load." * compensator_loads[i])
#         _ODSS.Loads.kW(P_comp)
#         _ODSS.Loads.kvar(Q_comp)

#         println("Updated Compensator Load for Phase $(phases[i]): P = $P_comp kW, Q = $Q_comp kvar")
#     end
# end

### not fro main network n
_ODSS.Circuit.SetActiveElement("Load.70")
currents_load4 = _ODSS.CktElement.CurrentsMagAng()
magnitude_2c1 = currents_load4[1]
angle_2c1 = currents_load4[2]
I2_c1 = round(magnitude_2c1 * cis(deg2rad(angle_2c1)),digits=2)
_ODSS.Circuit.SetActiveElement("Load.64")
currents_load5 = _ODSS.CktElement.CurrentsMagAng()
magnitude_2c2 = currents_load5[1]
angle_2c2 = currents_load5[2]
I2_c2 = round(magnitude_2c2 * cis(deg2rad(angle_2c2)),digits=2)
I2_c = I2_c1 + I2_c2

#########################################


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


# Plotting Load Currents (I2)
println("Plotting Load Currents (I2)")
phasorcosine(abs(I2_before[1]), angle(I2_before[1]), ylabel=L"$i$", maglabel=L"$\hat{I2}_a$", 
    labelrsep=0.5, figsize=(7*a, 2.5*a), color="blue", linestyle="-")
phasorcosine(abs(I2_before[2]), angle(I2_before[2]), ylabel=L"$i$", maglabel=L"$\hat{I2}_b$ ", 
    labelrsep=0.5, color="green", linestyle="-", add=true)
phasorcosine(abs(I2_before[3]), angle(I2_before[3]), ylabel=L"$i$", maglabel=L"$\hat{I2}_c$ ", 
    labelrsep=0.5, color="red", linestyle="-", add=true)
phasorcosine(abs(I2_n_before), angle(I2_n_before), ylabel=L"$i$", maglabel=L"$\hat{I2}_n$ ", 
    labelrsep=0.5, color="black", linestyle="--", add=true)
# Show the plot for Load Currents
gcf()
show()

println("Plotting Compensating Currents (I3)")
phasorcosine(abs(I3[1]), angle(I3[1]), ylabel=L"$i$", maglabel=L"$\hat{I3}_a$", 
    labelrsep=0.5, figsize=(7*a, 2.5*a), color="blue", linestyle="-")
phasorcosine(abs(I3[2]), angle(I3[2]), ylabel=L"$i$", maglabel=L"$\hat{I3}_b$", 
    labelrsep=0.5, color="green", linestyle="-", add=true)
phasorcosine(abs(I3[3]), angle(I3[3]), ylabel=L"$i$", maglabel=L"$\hat{I3}_c$ ", 
    labelrsep=0.5, color="red", linestyle="-", add=true)
phasorcosine(abs(I3_n), angle(I3_n), ylabel=L"$i$", maglabel=L"$\hat{I3}_c$ ", 
    labelrsep=0.5, color="black", linestyle="--", add=true)
# Show the plot for Compensating Currents
gcf()
show()

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

println("Plotting Substation Voltage before Compensation")
phasorcosine(abs(voltage_sub_bef[1]), angle(voltage_sub_bef[1]), ylabel=L"$i$", maglabel=L"$\hat{V1}_a$", 
    labelrsep=0.5, figsize=(7*a, 2.5*a), color="blue", linestyle="--")
phasorcosine(abs(voltage_sub_bef[2]), angle(voltage_sub_bef[2]), ylabel=L"$i$", maglabel=L"$\hat{V1}_b$", 
    labelrsep=0.5, color="green", linestyle="--", add=true)
phasorcosine(abs(voltage_sub_bef[3]), angle(voltage_sub_bef[3]), ylabel=L"$i$", maglabel=L"$\hat{V1}_c$ ", 
    labelrsep=0.5, color="red", linestyle="--", add=true)
phasorcosine(abs(voltage_sub_bef[4]), angle(voltage_sub_bef[4]), ylabel=L"$i$", maglabel=L"$\hat{V1}_n$ ", 
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


########### function to plot the unbalance factor

## TO ARRANGE THE LOADS BASED ON THE DISTANCE FROM THE SUBSTATION and THEIR PHASES
## Get load voltages
function get_load_voltages()
    load_voltages = Dict{String, Vector{Float64}}()
    for load_name in _ODSS.Loads.AllNames()
        _ODSS.Circuit.SetActiveElement("Load." * load_name)
        voltages = _ODSS.CktElement.VoltagesMagAng()
        load_voltages[load_name] = voltages[1:2:end]  # Collecting all phase voltages
    end
    return load_voltages
end
 load_voltages = get_load_voltages()



function get_bus_components(;comp_types=["line", "transformer", "load", "gen"])::Dict{String,Any}
    bus_components = Dict{String,Any}()
    for bus_name in _ODSS.Circuit.AllBusNames()
        _ODSS.Circuit.SetActiveBus(bus_name)

        pde_elements = _ODSS.Bus.AllPDEatBus()
        loads = _ODSS.Bus.LoadList()

        bus_components[bus_name] = Dict(comp_type => Any[] for comp_type in comp_types)
        for pde_element in pde_elements
            if startswith(pde_element, "Line")
                append!(bus_components[bus_name]["line"], split(pde_element,".")[2:end])
            elseif startswith(pde_element, "Transformer")
                append!(bus_components[bus_name]["transformer"], split(pde_element,".")[2:end])
            end
        end
        for load in loads
            append!(bus_components[bus_name]["load"], split(load,".")[2:end])
        end
    end
    return bus_components
end
 
bus_components = get_bus_components()

# Get line data
function get_line_data()
    line_data = Dict()
    for line_name in OpenDSSDirect.Lines.AllNames()
        OpenDSSDirect.Lines.Name(line_name)
        line_data[line_name] = Dict("length" => OpenDSSDirect.Lines.Length(),
                                    "bus1" => OpenDSSDirect.Lines.Bus1(),
                                    "bus2" => OpenDSSDirect.Lines.Bus2())
    end
    return line_data
end
line_data = get_line_data()

function get_bus_voltage_snap()
    bus_dict = Dict()
    phase_mapping = Dict(1 => "a", 2 => "b", 3 => "c", 4 => "n")

    for bus_name in OpenDSSDirect.Circuit.AllBusNames()
        bus_dict[bus_name] = Dict()
        OpenDSSDirect.Circuit.SetActiveBus(bus_name)

        bus_phases = OpenDSSDirect.Bus.Nodes()  # Get the connected phases of the bus
        bus_voltages = OpenDSSDirect.Bus.PuVoltage()  # Get per-unit voltages of the bus

        # Check if neutral (phase 4) is present
        has_neutral = 4 ∈ bus_phases
        if has_neutral
            neutral_voltage = bus_voltages[length(bus_voltages)] 
        else
            neutral_voltage = 0+0im
        end

        # Remove phase 4 from the list of bus phases for processing
        bus_phases = filter(x -> x != 4, bus_phases)

        for (i, phase) in enumerate(bus_phases)
            voltage_phase = bus_voltages[i]
            # Calculate voltage magnitude and angle, with or without neutral reference
            if has_neutral
                bus_dict[bus_name]["vm" * phase_mapping[phase]] = abs(voltage_phase - neutral_voltage)
                bus_dict[bus_name]["va" * phase_mapping[phase]] = angle(voltage_phase - neutral_voltage) * 180 / pi
            else
                bus_dict[bus_name]["vm" * phase_mapping[phase]] = abs(voltage_phase)
                bus_dict[bus_name]["va" * phase_mapping[phase]] = angle(voltage_phase) * 180 / pi
            end
        end
    end

    return bus_dict
end
# Measure system voltages
bus_dict = get_bus_voltage_snap()  # Get the voltages of all phases of all buses


## include the load, their distance and the connected bus
function calculate_load_distances_from_bus(starting_bus, bus_components, line_data)
    load_distances = Dict{String, Tuple{Float64, String}}()  # Stores distance and connected bus
    visited_buses = Set{String}() 
    function traverse_bus(bus, current_distance)
        if bus in visited_buses
            return
        end    
        # println("Traversing Bus: $bus, Distance from $starting_bus: $current_distance km")
        push!(visited_buses, bus)
        # Check for loads at the current bus
        if "load" in keys(bus_components[bus])
            for load in bus_components[bus]["load"]
                println("Found Load: $load at Bus: $bus, Distance: $current_distance km")
                load_distances[load] = (current_distance, bus)  # Store distance and bus name
            end
        end
        # Traverse connected lines
        if "line" in keys(bus_components[bus])
            for line_name in bus_components[bus]["line"]
                if line_name in keys(line_data)
                    line_info = line_data[line_name]
                    next_bus = nothing
                    if startswith(line_info["bus1"], bus)
                        next_bus = split(line_info["bus2"], ".")[1]
                    elseif startswith(line_info["bus2"], bus)
                        next_bus = split(line_info["bus1"], ".")[1]
                    end
                    if next_bus != nothing && !(next_bus in visited_buses)
                        traverse_bus(next_bus, current_distance + line_info["length"])
                    end
                end
            end
        end
    end
    traverse_bus(starting_bus, 0.0)
    return load_distances
end

load_distances = calculate_load_distances_from_bus("6687", bus_components, line_data)
#plot the load distances
distances = [v[1] for v in values(load_distances)]
load_names = collect(keys(load_distances))
# Sort by distance
sorted_indices = sortperm(distances)
sorted_distances = distances[sorted_indices]
sorted_load_names = load_names[sorted_indices]
# Scatter plot
scatter(
    sorted_distances,
    sorted_load_names,
    xlabel = "Distance (km)",
    ylabel = "Load Name",
    title = "Scatter Plot of Load Names vs Distance",
    legend = false,
    color = :blue,
    xticks = 0:0.05:0.4,
    markersize = 2
)
#################

using Plots


function plot_voltage_magnitudes_perphase_with_labels(load_distances, bus_dict)
    # Define phases and initialize plots
    phases = ["a", "b", "c"]
    for phase in phases
        distances = Float64[]
        voltages = Float64[]
        labels = String[]
        
        for (load, (distance, bus)) in load_distances
            phase_key = "vm" * phase  # Key for phase voltage
            if haskey(bus_dict[bus], phase_key)
                push!(distances, distance)
                push!(voltages, bus_dict[bus][phase_key])
                push!(labels, string(load))  # Use load or bus name as label
            end
        end
        # Sort data by distance for better visualization
        sorted_indices = sortperm(distances)
        distances = distances[sorted_indices]
        voltages = voltages[sorted_indices]
        labels = labels[sorted_indices] 
        # Plot voltage magnitudes
        p = plot(distances, voltages, label="Phase $phase", xlabel="Distance from Substation (km)",
            ylabel="Voltage Magnitude (p.u.)", title="Loads voltage Magnitude - Phase $phase", linewidth=2)

        # Add labels to points
        for (x, y, label) in zip(distances, voltages, labels)
            annotate!(x, y, text(label, :red, 8))
        end
        # Save or display the plot
        savefig(p, "loads voltage_magnitude_phase_$phase.png")
        display(p)
    end
end
# Call the function
plot_voltage_magnitudes_perphase_with_labels(load_distances, bus_dict)



function plot_voltage_magnitudes_separate_phases_with_labels(load_distances, bus_dict_before, bus_dict_after)
    """
    Plot voltage magnitudes before and after adjustments for each phase in separate plots,
    with load names labeled on the plot.

    Args:
        load_distances: Dictionary of loads with distances and connected buses.
        bus_dict_before: Dictionary of voltage data for buses before adjustments.
        bus_dict_after: Dictionary of voltage data for buses after adjustments.
    """
    # Define phases
    phases = ["a", "b", "c"]
    
    for phase in phases
        distances = Float64[]
        voltages_before = Float64[]
        voltages_after = Float64[]
        labels = String[]
        
        for (load, (distance, bus)) in load_distances
            phase_key = "vm" * phase  # Key for phase voltage
            
            # Check if the phase exists in both before and after data
            if haskey(bus_dict_before[bus], phase_key) && haskey(bus_dict_after[bus], phase_key)
                push!(distances, distance)
                push!(voltages_before, bus_dict_before[bus][phase_key])
                push!(voltages_after, bus_dict_after[bus][phase_key])
                push!(labels, string(load))  # Store load name for annotation
            end
        end
        
        # Sort data by distance for better visualization
        sorted_indices = sortperm(distances)
        distances = distances[sorted_indices]
        voltages_before = voltages_before[sorted_indices]
        voltages_after = voltages_after[sorted_indices]
        labels = labels[sorted_indices]
        
        # Create a separate plot for each phase
        p = plot(distances, voltages_before, label="Before Adjustment", linewidth=2, linestyle=:solid,
                 xlabel="Distance from Substation (km)", ylabel="Voltage Magnitude (p.u.)",
                 title="Voltage Magnitudes - Phase $phase", legend=:topright)
        plot!(p, distances, voltages_after, label="After Adjustment", linewidth=2, linestyle=:dash)

        # Add labels to the points
        for (x, y, label) in zip(distances, voltages_before, labels)
            annotate!(p, x, y, text(label, :blue, 8))  # Label for "before"
        end
        for (x, y, label) in zip(distances, voltages_after, labels)
            annotate!(p, x, y, text(label, :red, 8))  # Label for "after"
        end

        # Save and display the plot
        savefig(p, "voltage_magnitudes_phase_$phase.png")
        display(p)
    end
end

# Example call to the function
plot_voltage_magnitudes_separate_phases_with_labels(load_distances, bus_dict1, bus_dict2)



##

# Filter buses with exactly six elements to calculate the unbalance
bus_dict_3phase1 = Dict(k => v for (k, v) in bus_dict1 if length(keys(v)) == 6)
bus_dict_3phase2 = Dict(k => v for (k, v) in bus_dict2 if length(keys(v)) == 6)

function measure_voltage_unbalance(filtered_bus_dict, calculate_symmetrical_components)
    unbalance_dict = Dict()
    for (bus, voltages) in filtered_bus_dict
        # Extract phase voltages
        vma = voltages["vma"]
        vaa = deg2rad(voltages["vaa"])  # Convert angle to radians
        vmb = voltages["vmb"]
        vab = deg2rad(voltages["vab"])
        vmc = voltages["vmc"]
        vac = deg2rad(voltages["vac"])  
        # Create voltage phasors
        Va = vma * cis(vaa)
        Vb = vmb * cis(vab)
        Vc = vmc * cis(vac) 
        # Calculate symmetrical components
        symmetrical_components = calculate_symmetrical_components([Va, Vb, Vc])
        V0, V1, V2 = symmetrical_components  # Zero, positive, and negative sequence voltages
        # Compute voltage unbalance factor
        # unbalance_factor = sqrt(abs(V0)^2 + abs(V2)^2) / abs(V1)*100
        unbalance_factor = (abs(V2) / abs(V1)) * 100
        # Store the unbalance factor
        unbalance_dict[bus] = unbalance_factor
    end
    return unbalance_dict
end
unbalance_dict1 = measure_voltage_unbalance(bus_dict_3phase1, calculate_symmetrical_components)
unbalance_dict2 = measure_voltage_unbalance(bus_dict_3phase2, calculate_symmetrical_components)
# Find the bus with the highest unbalance
most_unbalanced_bus1 = argmax(unbalance_dict1)
most_unbalanced_bus2 = argmax(unbalance_dict2)

bus_numbers1 = collect(keys(unbalance_dict1))
unbalance_factors1 = collect(values(unbalance_dict1))
bus_numbers2 = collect(keys(unbalance_dict2))
unbalance_factors2 = collect(values(unbalance_dict2))
# Plotting
bar(bus_numbers1, unbalance_factors1,
    xlabel = "Bus Number",
    ylabel = "Voltage Unbalance Factor",
    title = "Voltage Unbalance Factor by Bus Number_before",
    bar_width = 0.6,
    legend = false,
    color = :orange,
    xticks = (1:length(bus_numbers1), bus_numbers1),
    xrotation=45,
    yticks = 0:0.5:4
)

bar(bus_numbers2, unbalance_factors2,
    xlabel = "Bus Number",
    ylabel = "Voltage Unbalance Factor",
    title = "Voltage Unbalance Factor by Bus Number_after",
    bar_width = 0.6,
    legend = false,
    color = :skyblue,
    xticks = (1:length(bus_numbers2), bus_numbers2),
    xrotation=45,
    yticks = 0:0.5:6.5
)

# ##Combine data into a grouped bar plot
bus_numbers = bus_numbers1  # Assumes both datasets share the same bus numbers
# Offset for clearer distinction
offset = 0.2
# Bar positions
x_positions1 = collect(1:length(bus_numbers)) .- offset
x_positions2 = collect(1:length(bus_numbers)) .+ offset
# Plot
bar(x_positions1, unbalance_factors1,
    width = 0.4,
    label = "Before",
    color = :orange,
    edgecolor = :black,
    xticks = (1:length(bus_numbers1), bus_numbers1),
    xrotation = 45,
)
bar!(x_positions2, unbalance_factors2,
    width = 0.2,
    label = "After",
    color = :skyblue,  # Set bar fill color to white for outline visibility
    linecolor = :orange,
    linestyle = :dash,
    linewidth = 2
)

# Customization
xlabel!("Bus Number")
ylabel!("Voltage Unbalance Factor")
title!("Voltage Unbalance Factor by Bus Number")


# Display plot


#### for single phase buses
filtered_bus_dict2 = Dict(k => v for (k, v) in bus_dict if length(keys(v)) == 2)
v_nominal = 1.0
# Calculate unbalance for each bus
unbalance_dict = Dict{String, Float64}()
for (bus, values) in filtered_bus_dict2
    for (key, value) in values
        if startswith(key, "vm")  # Only process voltage magnitudes
            unbalance = abs(value - v_nominal) / v_nominal * 100
            unbalance_dict[bus] = unbalance
        end
    end
end
# Extract data for plotting
buses = collect(keys(unbalance_dict))
unbalances = collect(values(unbalance_dict))
# Sort by unbalance
sorted_indices = sortperm(unbalances)
sorted_buses = buses[sorted_indices]
sorted_unbalances = unbalances[sorted_indices]
# Plot
bar(sorted_buses, sorted_unbalances,
    xlabel = "Bus numbers",
    ylabel = "Unbalance(%), v_phase-v_nom/v_nom",
    title = "Voltage Unbalance for Single-Phase Buses",
    bar_width = 0.6,
    legend = false,
    color = :orange,
    xtickfont = font(4),  # Set smaller font size for bus numbers
    xticks = (1:length(sorted_buses), sorted_buses),
    xrotation=45,  # Show all bus numbers
    yticks = 0:2:20
)

#########

######################





function get_currents()
    _ODSS.Circuit.SetActiveElement("Line.line_75588076_6694_6625")
    currents = _ODSS.CktElement.CurrentsMagAng()
    I_a = round(currents[1,5] * cis(deg2rad(currents[2,5])), digits=2)
    I_b = round(currents[1,6] * cis(deg2rad(currents[2,6])), digits=2)
    I_c = round(currents[1,7] * cis(deg2rad(currents[2,7])), digits=2)
    I_n = round(currents[1,8] * cis(deg2rad(currents[2,8])), digits=2)
    return [I_a, I_b, I_c, I_n]
end

I2_bef= get_currents()
I2_a_before, I2_b_before, I2_c_before, I2_n_before= I2_bef
I2_before = I2_bef[1:3]




### just for testing
all_elements = OpenDSSDirect.Circuit.AllElementNames()## returns all th elements in the network
element_name = "Load.1"  # Load.1 ## this way we see this element is connected to what
_ODSS.Circuit.SetActiveElement(element_name)
connected_buses = _ODSS.CktElement.BusNames()
_ODSS.Circuit.SetActiveBus("6687")
voltage = _ODSS.Bus.PuVoltage()
println("Bus Voltage (pu): ", voltage)
_ODSS.Circuit.SetActiveElement("Load.68")
load68_powers = _ODSS.CktElement.Powers()
println("Load 68 powers: ", load68_powers)
_ODSS.Circuit.SetActiveElement("Transformer.tr1")
transformer_currents = _ODSS.CktElement.CurrentsMagAng()
println("Transformer currents: ", transformer_currents)