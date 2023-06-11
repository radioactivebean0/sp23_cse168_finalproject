# To-dos

- More reading

## Project

- Fix CLI
- KD tree support three operations:
	- Keys: 3 doubles (x, y, z) world coordinates, no value
	- One thread (but look into multiple): construction
	- Gather contributions of all photons near visible points
		- Count number of nodes in KD tree within a sphere of radius R of a given point, don't care
		- about coordinates
	- Clear

## Goal

### June 10th

Must
- First pass of algorithm done; trace ray from camera into the scene
- Integrating KD-tree lib
- One different scene for debugging (at least point light with cbox water, possibly with a different
  angle)

Optional
- Placing photons in the scene
- One different model of cup/glass for caustics
- A scene of pouring water into ice
	- <https://benedikt-bitterli.me/resources/images/glass-of-water.png>
