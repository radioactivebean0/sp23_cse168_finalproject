import math
y = 0.0
z = -2.0
r = 1
numPoints = 26
sphere_r = r * math.sin(math.pi/numPoints)
points = []
for index in range(numPoints):
    points.append((r*math.cos((index*2*math.pi)/numPoints), y , r*math.sin((index*2*math.pi)/numPoints)-2.0))

for pt in points:
    print(f'<shape type="sphere"> <point name="center" x="{pt[0]}" y="{pt[1]}" z="{pt[2]}"/> <float name="radius" value="{sphere_r}"/></shape>')
