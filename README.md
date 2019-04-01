# Open-Channel-Flow

## installation
```
virtualenv env
source env/bin/activate

```



```
pip install -r requirements.txt

```




```
from geometry import ChannelSection, FlowArea

coordinates = (
    (0, 0),
    (0,-10),
    (60, -10),
    (60, 0),
)

A = ChannelSection(coordinates)
A.plot_geometry()

Af = FlowArea(A, 2.5, 2500)    
Af.plot_geometry()

print("Centroid :", Af.get_centroid())

print("Specific Force :", Af.specific_force())

print("velocity : ", Af.velocity())

print("Top Width: ", Af.free_surface_length())
print("Hydraulic Radius : ", Af.hydraulic_radius())
print("Hydraulic Depth : ", Af.hydraulic_depth())
print("Hydraulic Jump y2 : ", Af.hydraulic_jump_y2())

print("Area cross-section flow:",  round(Af.area(), 2))
print("Froud number :", Af.froude_number())
Af.specific_energy_plot()

```
![alt text](https://raw.githubusercontent.com/enkaynitin/Open-Channel-Flow/master/Channel%20Section.png)
![alt text](https://raw.githubusercontent.com/enkaynitin/Open-Channel-Flow/master/Channel%20Section%20with%20Flow.png)
![alt text](https://raw.githubusercontent.com/enkaynitin/Open-Channel-Flow/master/Specific%20Energy.png)