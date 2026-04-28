import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.patches import Polygon
import yaml

# Define parameters from yaml file
with open("src/PF/ParticleFilterParams.yaml", "r") as file:
    config = yaml.safe_load(file)

# Define problem definition yaml file
with open("input/problem.yaml", "r") as file:
    problemConfig = yaml.safe_load(file)

maxX = problemConfig["x_max"]
maxY = problemConfig["y_max"]
fov = np.pi/config["sensor"]["fovFractionOfPi"]
sensorRange = config["sensor"]["range"]
gridSizeX = config["Filter"]["gridSizeX"]
gridSizeY = config["Filter"]["gridSizeY"]
observedDecayRate = config["Filter"]["observedDecayRate"]
maxObservedScore = config["Filter"]["maxObservedScore"]

# Define obstacles - make sure these are the same as those used in main.cpp
obstacles = []
for obstacleConfig in problemConfig["obstacles"]:
    vertices = []
    for vertexConfig in obstacleConfig["vertices_ccw"]:
        if isinstance(vertexConfig, dict):
            vertices.append([vertexConfig["x"], vertexConfig["y"]])
        else:
            vertices.append([vertexConfig[0], vertexConfig[1]])
    obstacles.append(np.array(vertices, dtype=float))

# Load the most recent output csv file
data = np.loadtxt(config["misc"]["outputFile"], delimiter=",")

# Figure setup
fig, ax = plt.subplots(figsize=(7,7))
ax.set_xlim(0, maxX)
ax.set_ylim(0, maxY)
ax.set_aspect("equal")
ax.set_xlabel("X Position")
ax.set_ylabel("Y Position")

# Create a heatmap for visualizing the observed grid
observedGrid = np.zeros((gridSizeX, gridSizeY))
heatMap = ax.imshow(observedGrid.T, extent=[0, maxX, 0, maxY], origin="lower", cmap="coolwarm", vmin=0, vmax=maxObservedScore, alpha=0.5, interpolation="nearest")

# Create a colorbar
colorBar = plt.colorbar(heatMap, ax=ax)
colorBar.set_label("Observed Score")

# Create obstacles in figure
for obstacle in obstacles:
    shape = Polygon(obstacle, closed=True, facecolor="gray", edgecolor="black", alpha=1)
    ax.add_patch(shape)

# Plot the particles, vehicle, target, and filter estimate over time
particlesPlot, = ax.plot([], [], ".", color="blue", markersize=2, label="Particles")
vehiclePlot, = ax.plot([], [], "o", color="lime", markersize=10, label="Vehicle")
targetPlot, = ax.plot([], [], "x", color="red", markersize=10, label="Target")
estimatePlot, = ax.plot([], [], "*", color="lime", markersize=10, label="Estimate")
densityEstimatePlot, = ax.plot([], [], "s", color="cyan", markersize=10, label="Density Estimate")  

# Plot the field of view path over time
fovPatch = Polygon([[0,0]], closed=True, color="lime", alpha=0.2, label="Sensor FOV")
ax.add_patch(fovPatch)

# Add a legend to the plot
ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.05), ncol=3)
plt.tight_layout()

# Checks if points are in counter clockwise order
def counterClockwise(p1, p2, p3):
    return (p3[1] - p1[1])*(p2[0] - p1[0]) > (p2[1] - p1[1])*(p3[0] - p1[0])

# Checks if the two segments intersect
def segmentsIntersect(p1, p2, p3, p4):
    return (counterClockwise(p1, p3, p4) != counterClockwise(p2, p3, p4)) and (counterClockwise(p1, p2, p3) != counterClockwise(p1, p2, p4))

# Check whether the vehicle can observe the provided grid cell location (if there is an obstacle in the way)
def lineOfSight(vehicle, target):
    for obstacle in obstacles:
        for i in range(len(obstacle)):
            obstacleVertex1 = obstacle[i]
            obstacleVertex2 = obstacle[(i+1)%len(obstacle)]
            if segmentsIntersect(vehicle, target, obstacleVertex1, obstacleVertex2):
                return False
    return True

# Calculate the fov shape of the sensor
def fovShape(vehicle, heading):
    headingAngle = np.arctan2(heading[1], heading[0])
    fovVertexAngles = np.linspace(headingAngle - fov/2, headingAngle + fov/2, 50)

    fovShapeVertices = [vehicle]
    for angle in fovVertexAngles:
        fovShapeVertices.append([vehicle[0] + sensorRange*np.cos(angle), vehicle[1] + sensorRange*np.sin(angle)])

    return np.array(fovShapeVertices)

# Check if a grid cell is observed
def observed(cellPosition, vehicle, heading):
    dx = cellPosition[0] - vehicle[0]
    dy = cellPosition[1] - vehicle[1]

    distanceSquared = dx*dx + dy*dy
    if distanceSquared > sensorRange**2: return False

    dot = dx*heading[0] + dy*heading[1]

    if dot <= np.sqrt(distanceSquared)*np.cos(fov/2.0):
        return False

    if not lineOfSight(vehicle, cellPosition):
        return False

    return True

# Update the frame of the animation
def updateFrame(frameIndex):
    global observedGrid

    # Get the positions of all objects of interest from the current row of the csv file
    currentPositions = data[frameIndex]
    vehicle = currentPositions[0:2]
    heading = currentPositions[2:4]
    estimate = currentPositions[4:6]
    target = currentPositions[6:8]
    estimate_density = currentPositions[8:10]
    particles = currentPositions[10:].reshape(-1,2)

    # Update the plots with the current vehicle, estimate, target, and particle positions
    vehiclePlot.set_data([vehicle[0]], [vehicle[1]])
    targetPlot.set_data([target[0]], [target[1]])
    estimatePlot.set_data([estimate[0]], [estimate[1]])
    densityEstimatePlot.set_data([estimate_density[0]], [estimate_density[1]])
    particlesPlot.set_data(particles[:,0], particles[:,1])

    # Update the sensor fov shape
    fovPatch.set_xy(fovShape(vehicle, heading))

    # Update the observed grid by decaying the values and then increase the score of the currently observed grid cells
    observedGrid *= observedDecayRate

    for i in range(gridSizeX):
        for j in range(gridSizeY):
            gridCell = np.array([(i + 0.5)*maxX/gridSizeX, (j + 0.5)*maxY/gridSizeY])

            if observed(gridCell, vehicle, heading):
                observedGrid[i, j] = min(observedGrid[i, j] + 1.0, maxObservedScore)

    # Redefine the heatmap based on the current observed grid
    heatMap.set_array(observedGrid.T)

    # Dynamically update the plot title
    ax.set_title(f"Particle Filter Estimate at Frame {frameIndex}")

    return (particlesPlot, vehiclePlot, targetPlot, estimatePlot, densityEstimatePlot, fovPatch, heatMap)

# Create animation
animation = FuncAnimation(fig, updateFrame, frames=data.shape[0], interval=5, repeat=False)

# Display the animation
plt.show()

# Save animation video as a file
# animation.save("particle_filter.mp4", fps=5)