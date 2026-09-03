“Simple” plots for PDAC modeling paper

Note: All these figures should be generated using the pdacagviz module in article sizing

1. Which classes of tumor cells grow out under each condition/ each tissue setting
- scatter plot
- x-axis as either inital number of cancer cells or final number of cancer cells
- y-axis as percentage change (final-inital)/inital*100 of cancer cells
- each dot is a simulation run
- dot should be colored by simulation type
- three columns, each for different cancer cell type (type1, type2, type1+type2)
- x and y axis should have same min/max limits between the three columns for comparison
- a fourth scatter with x-axis as percentage change of type1 and y-axis as percentage change of type1+type2

2. How many CD8 attacks happen in total
Good to know how many of these were against tumor expressing just class I vs class I and class II
Check this for each simulation set
May want to check between all sets, normalizing for the cell count
- let's first explore the derived data to see if extracting this interaction is possible

3. htan_geometries will be a challenge to design visualizations for but I think if we are smart it can be really simple/clean

4. How many Tregs induced (by apCAF contact) and does having more Tregs by this mechanism associate with less T cell kills, different tumor growth, etc
	•		•	This is also where we would compare the different model mechanisms, maybe we can submit those other jobs?