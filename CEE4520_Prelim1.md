Name: Allison Tran (ant42)
CEE 4520 Prelim 1

# Prelim 1 CEE 4520
Open book, open internet. No discussion boards, no posting issues, no conversation with anyone other than Monroe, Matt, or Jillian about anything on this exam! You may send questions to Monroe, Matt, and Jillian.

Note that total points for the exam = 106

# Multiple choice (4 points each)

**Bold** the correct answer

 1) How far away from the location where the chemical dose is set is the flocculator in an AguaClara plant (pick the answer that is the best approximation)?

    **a) 1 m**

    b) 5 m

    c) 10 m

    d) 20 m

2) Which equation can be used when there is a loss (or increase) of mechanical energy?

    **a) Energy Equation**

    b) Bernoulli Equation

3) What happens to the maximum energy dissipation rate in a plane jet if velocity is held constant and the width of the jet is increased?

    a) The maximum energy dissipation rate increases

    b) The maximum energy dissipation rate remains constant

    **c) The maximum energy dissipation rate decreases**

4) Which type of head loss increases linearly with flow rate?

    **a) Laminar flow major loss**

    b) Laminar flow minor loss

    c) Turbulent flow major loss

    d) Turbulent flow minor loss

# Short Answer (4 points each)

1) What is the dimensionless *vena contracta* coefficient a ratio of?

**The vena contracta coefficient is the ratio of the area of the jet exiting the orifice to the area of the orifice.**


2) Identify a location in an AguaClara entrance tank where the Bernoulli equation can be used to describe the flow. Explain WHY the Bernoulli equation can be used.

**The Bernoulli equation is applicable between the exit of the jet diffuser and the exit of the jet reverser because there is no flow expansion yet.**

3) What happens to the coagulant dose (assuming the operator doesn't change anything) in an AguaClara plant when the flow rate through the plant changes from the plant maximum design flow to 50% of the maximum design flow?

**The coagulant dose should also decrease by 50% because the relationship between the coagulant dose and the design flow rate should be linear.**


4) What is Monroe's hypothesis for why primary particles such as clay and pathogens can't attach to flocs in a flocculator?

**Monroe's hypothesis for why primary particles (such as clay and pathogens) can't attach to flocs in a flocculator is because of the boundary layer of fluid around the floc.**

5)  Give one example each of a minor loss and a major loss in AguaClara plants

**One source of minor loss in AguaClara plants are baffles, which are obstructions in the flow of the water to switch the direction of the flow by 180 degrees. Baffles are also a source of major loss in AguaClara plants, which can come from changes in flocculation elevation.**



# Design Challenges
Document your work. Show any equations that you use in latex. Then do the calculation in python. Simplify the units to something that is easy to understand (SI metric system please!). Print your answer in python and copy that answer and paste it into markdown so your answer is shown without needing to run the code.

Here are some import statements that will likely be helpful. You may, of course, import other packages if needed!

```Python
from aguaclara.core.units import unit_registry as u
from aguaclara.core import utility as ut
from aguaclara.core import physchem as pc
import aguaclara.design.floc as floc
from aguaclara.research import floc_model as fm
import matplotlib.pyplot as plt
import numpy as np
```

1) Water at $15^\circ C$ flows through a straight, smooth tube that is 1/8th inch in inner diameter, 2 m long, has a head loss of 20 cm, and has total minor loss coefficients for entrance and exit of 4.
a) (5 points) What is the flow rate?
b) (5 points) Is this flow laminar or turbulent?

```Python
Water_temp = (15*u.degC).to(u.degK)
Inner_diameter = (0.125*u.inch).to(u.m)
Length = 2*u.m
Head_loss = (20*u.cm)
Minor_loss_coefficient = 4
Nu = pc.viscosity_kinematic(Water_temp) #solving for the kinematic viscosity
PipeRough = 0 #roughness of a smooth pipe is 0

#solving for the flow rate
Flow_rate = pc.flow_pipe(Inner_diameter, Head_loss, Length, Nu, PipeRough, Minor_loss_coefficient)
Flow_rate = Flow_rate.to(u.L/u.sec)

print('The flow rate is',ut.round_sf(Flow_rate,3))

#solving for the Reynolds number to determine laminar or turbulent
Rey_Num = pc.re_pipe(Flow_rate, Inner_diameter, Nu)
Rey_Num

print('The flow is laminar because the Reynolds number for laminar flow is less than 2100, and the Reynolds number of this flow is',ut.round_sf(Rey_Num,3))
```
**The flow rate is 0.00201 liters per second.
The flow is laminar because the Reynolds number for laminar flow is less than 2100, and the Reynolds number of this flow is 709.**

2) (10 points) A mechanical hydraulic mix unit for a 3 mgd (million gallon per day) water treatment plant has a residence time of 30 seconds and a $G_{CS}$ of 1000 Hz. Estimate the cost of this energy per year.  You may neglect the fact that the designer had to round up to next available motor size. You may assume that the rapid mix unit is 80% efficient at converting electricity into fluid shear. The temperature ranges from $0 \circ C$ to $30 \circ C$. The unit was designed to deliver the target fluid deformation for the worst case temperature. The cost of electricity is 0.15 USD/(kW hr).

Equation for energy dissipation rate: $G_{CS} = \sqrt{\frac{\bar \varepsilon}{\nu}}$

```python
Q = 3000000*(u.gallon/u.day)
Residence_time = 30*u.sec
G_CS = 1000*u.Hz
Temp_low = 0*u.degC
Temp_high = 30*u.degC

#solving for dynamic viscosity at minimum and maximum temperature
Nu_low = pc.viscosity_dynamic(Temp_low)
Nu_high = pc.viscosity_dynamic(Temp_high)

#solving for the energy rate at low temperature
energy_rate_low = (G_CS**2)*Nu_low
energy_rate_low = energy_rate_low.to_base_units()
energy_rate_low

#converting to fluid shear from the energy rate
fluid_to_energy_rate_low = energy_rate_low*(Q.to(u.m**3/u.sec))*Residence_time/0.8 #0.8 is efficiency of conversion
fluid_to_energy_rate_low
fluid_to_energy_rate_low = fluid_to_energy_rate_low.to(u.joules/u.sec)
fluid_to_energy_rate_low

#calculating the cost of energy per year at the given fluid shear
price = 0.15*(u.dollars/(u.kilowatthour))
time = (1*u.year).to(u.sec)
cost_low = (price * fluid_to_energy_rate_low) * time
cost_low = cost_low.to_base_units()
cost_low

#solving for the energy rate at high temperature
energy_rate_high = (G_CS**2)*Nu_high
energy_rate_high = energy_rate_high.to_base_units()
energy_rate_high

fluid_to_energy_rate_high = energy_rate_high*(Q.to(u.m**3/u.sec))*Residence_time/0.8
fluid_to_energy_rate_high = fluid_to_energy_rate_high.to(u.joules/u.sec)
fluid_to_energy_rate_high

cost_high = (price*fluid_to_energy_rate_high) * time
cost_high = cost_high.to_base_units()
cost_high

print('The cost of energy would be',ut.round_sf(cost_low,6),'per year')
```
**The cost of energy per year would be $11,361 accommodating for a temperature of 0 degrees Celsius.**

3) A flocculator designed using the AguaClara code has a design flow of 20 L/s and a design temperature of $15^\circ C$.  You may use the flocculator code to check your answer, but do the calculation using an equation that you derive.

   a) (5 points) Calculate the Camp Stein velocity gradient, $G_{CS}$, for the flocculator based on the flow characteristics and geometry of one flow expansion.

   b) (5 points) Calculate the residence time for one flow expansion based on the flow characteristics and geometry of that expansion.

   c) (5 points) Calculate the $G_{CS}\theta$ for one flow expansion

   d) (5 points) Calculate the number of flow expansions required given the design goal of $G_{CS}\theta = 37,000 $.

Equation for Camp-Stein velocity gradient: $G_{CS} = \sqrt{\frac{\bar \varepsilon}{\nu}}$

Equation for residence time: $\theta = \frac{gh_{\rm{L}}}{\bar\varepsilon} \qquad\qquad$


```Python
#given default parameters
Q_design = 20*u.L/u.s
temp = 15*u.degC
max_L = 6*u.m
Gt = 37000
HL = 40*u.cm
downstream_H = 2*u.m
ent_tank_L = 1.5*u.m
max_W = 42*u.inch
drain_t = 30*u.min
myF = floc.Flocculator(Q=Q_design, temp=15 * u.degC)
gravity = 9.8*(u.m/(u.sec**2))
myF.design

#solving for Camp Stein velocity gradient
Nu_2 = pc.viscosity_kinematic(temp)
E_bar = (gravity*HL)/myF.retention_time
GCS = ((E_bar/Nu_2)**(1/2)).to_base_units()
print('The Camp Stein velocity gradient is',ut.round_sf(GCS,3))

#solving for residence time
Residence_time = (gravity*HL)/E_bar
print('The residence time is',ut.round_sf(Residence_time,3))

#solving for Gcs*theta
G_CS_theta = GCS*Residence_time
G_CS_theta
print('The dimensionless value of the Camp Stein velocity gradient for one flow expansion is',ut.round_sf(G_CS_theta,6))

#solving for head loss per baffle
import math
K = floc.Flocculator.BAFFLE_K
g = 9.8*(u.m/(u.sec**2))
A = myF.channel_W*(ut.round_sf(myF.baffle_S,2))
V = Q_design/A
HL_perbaffle= (K*(V**2)/(2*g))
HL_perbaffle = HL_perbaffle.to(u.mm)
Num_baffle = HL/HL_perbaffle
Num_baffle = Num_baffle.to_base_units()
Num_baffle

Num_expansions = Num_baffle*myF.expansion_n
Num_expansions

print('The number of flow expansions required would be',ut.round_sf(Num_expansions,3))
```
**The Camp Stein velocity gradient is 93.2/second.
The residence time is 397 seconds.
The dimensionless value of the Camp Stein velocity gradient for one flow expansion is 36987.
The number of flow expansions required would be 154.**

4) (15 points) Plot $G_{CS}\theta$ for the flocculator in problem 3 if it is built and then operated over a range of flows from 25% to 100% of the design flow. Start by deriving an equation for $G_{CS}\theta$ that accounts for the fact that the flocculator is built (geometry is constant) and the flow rate is varying.

```Python
Q_design = 20*u.L/u.s
Q_25 = Q_design*0.25
Q_values = np.arange((Q_25).magnitude, Q_design.magnitude, 1)
G_CS_theta_values = np.zeros([len(Q_values)])

temp = 15*u.degC
Nu_x = pc.viscosity_kinematic(temp)

for x in range(Q_values.size):
  Q_x = (Q_values[x])
  myFx = floc.Flocculator(Q=Q_x, temp=temp)
  E_bar_x = (gravity*(myFx.HL)/(myFx.retention_time))
  G_CS_x = (((E_bar_x)/Nu_x)**(1/2)).to_base_units()
  G_CS_theta_values[x] = ((G_CS_x)*(myFx.retention_time)).magnitude

plt.plot(Q_values, G_CS_theta_values)
plt.xlabel('Flow (L/s)')
plt.ylabel('GCS*Theta (dimensionless)')
plt.title('GCS*Theta vs Flow Rate')
plt.show()
```
**The plot of GCSTheta vs Flow Rate is a horizontal line because GCSTheta is independent of the flow rate, as seen in the calculations done for GCSTheta within the for-loop.**

5) One of our goals it to produce water that is dramatically cleaner than regulations require.

    a) (5 points) How far apart will clay particles (equivalent spherical diameter of $5 \mu m$ and density of $2650 \frac{kg}{m^3}$) be if their concentration is 0.01 mg/L?

    b) (5 points) What is the mass of a clay particle?

    c) (5 points) How many clay particles would there be in a 250 mL glass of water?

Equation for separation distance: $\bar \Lambda  = \frac{1}{n_P^{\frac{1}{3}}} = {\rlap{-} V_{\rm{Surround}}}^\frac{1}{3}$

```Python
import numpy as np
from aguaclara.core.units import unit_registry as u
from aguaclara.core import physchem as pc, utility as ut
u.enable_contexts('chem')

density = 2650*(u.kg/u.m**3)
conc = 0.01*(u.mg/u.L)
diam = 5*u.micrometer

lambda_clay = ((density/conc)*((np.pi*diam**3)/6))**(1/3)
lambda_clay = lambda_clay.to_base_units()
lambda_clay = lambda_clay.to(u.mm)

print('The separation distance of these two clay particles will be', ut.round_sf(lambda_clay,3))

volume_particle = (4/3)*np.pi*((diam/2)**3)
mass_particle = (density*volume_particle).to(u.mg)

print('The mass of a clay particle is',ut.round_sf(mass_particle.to(u.nanogram),3))

num_particles = ((((250*u.mL).to(u.L))*conc)/mass_particle)

print('There would be',ut.round_sf(num_particles,5),'in a 250 mL glass of water')
```
**The separation distance of these two clay particles will be 2.59 millimeters.
The mass of a clay particle is 0.173 nanograms.
There would be 14410 particles in a 250mL glass of water.**
