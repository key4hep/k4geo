### Small script showing how geometry can be saved to root files, 
### that can be easily viewed using e.g https://root.cern/js/latest/
### as described in https://hep-fcc.github.io/fcc-tutorials/full-detector-simulations/Visualization/Visualization.html#detector-geometry

chmod u+x scripts/utils/dd4hep2root

# IDEA vertex
./scripts/utils/dd4hep2root -c FCCee_IDEA_Empty_o1_v01.xml Vertex_IDEA_o1_v01.xml -o fccee_idea_vertex.root

# Complete IDEA
./scripts/utils/dd4hep2root -c FCCee_IDEA_o1_v01.xml -o fccee_idea.root
