$cfg="ecoh.yaml"
julia ecoh.jl $cfg
python ../plot_linear.py $cfg
python ../plot_2d.py $cfg
python plot_dyn.py $cfg