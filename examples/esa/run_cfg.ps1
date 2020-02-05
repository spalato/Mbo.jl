$cfg="basic_ia.yaml"
julia basic_ia.jl $cfg
python ../plot_linear.py $cfg
python ../plot_2d.py $cfg
#python plot_dyn.py $cfg