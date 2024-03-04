#
install(FILES fixed_boundary/fixed_boundary_ex1.ipynb DESTINATION examples/TokaMaker/fixed_boundary COMPONENT examples)
install(FILES fixed_boundary/fixed_boundary_ex2.ipynb DESTINATION examples/TokaMaker/fixed_boundary COMPONENT examples)
install(FILES fixed_boundary/gNT_example DESTINATION examples/TokaMaker/fixed_boundary COMPONENT examples)
#
install(FILES HBT/HBT_geom.json DESTINATION examples/TokaMaker/HBT COMPONENT examples)
install(FILES HBT/HBT_mesh_ex.ipynb DESTINATION examples/TokaMaker/HBT COMPONENT examples)
install(FILES HBT/HBT_eq_ex.ipynb DESTINATION examples/TokaMaker/HBT COMPONENT examples)
install(FILES HBT/HBT_mesh.h5 DESTINATION examples/TokaMaker/HBT COMPONENT examples)
#
install(FILES LTX/LTX_geom.json DESTINATION examples/TokaMaker/LTX COMPONENT examples)
install(FILES LTX/LTX_mesh_ex.ipynb DESTINATION examples/TokaMaker/LTX COMPONENT examples)
install(FILES LTX/LTX_eq_ex.ipynb DESTINATION examples/TokaMaker/LTX COMPONENT examples)
# install(FILES LTX/LTX_mesh.h5 DESTINATION examples/TokaMaker/LTX COMPONENT examples)
#
install(FILES ITER/ITER_mesh_ex.ipynb DESTINATION examples/TokaMaker/ITER COMPONENT examples)
install(FILES ITER/ITER_baseline_ex.ipynb DESTINATION examples/TokaMaker/ITER COMPONENT examples)
install(FILES ITER/ITER_geom.json DESTINATION examples/TokaMaker/ITER COMPONENT examples)
install(FILES ITER/ITER_mesh.h5 DESTINATION examples/TokaMaker/ITER COMPONENT examples)