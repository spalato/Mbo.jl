Set-PSDebug -Trace 1
$cfg=$args[0]
$root=$(Select-String -Path $cfg -Pattern '^rootname:\W*(.*)$').Matches.Groups[1].value
julia .\third_order_xc.jl $cfg
ls ${root}_r[nr].bin | %{julia .\spec_ft.jl $cfg $_ $_.Name.Replace("_r", "_s")}
ls ${root}_s[nr].bin | %{python .\plot_spec.py $cfg $_}
julia .\spec_abs.jl $cfg $root ${root}_sa.bin
#rm *_[sr][rn].bin
python .\plot_spec.py $cfg ${root}_sa.bin
#python .\beatmap.py $cfg ${root}_sa.bin