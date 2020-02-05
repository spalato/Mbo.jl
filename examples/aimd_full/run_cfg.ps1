$cfg=$args
if (!$cfg) {
    $cfg=$(ls *.yaml)
}
foreach ($c in $cfg) {
    $root=$(Select-String -Path $cfg -Pattern '^rootname:\W*(.*)$').Matches.Groups[1].value
    julia .\linear_xcf.jl $cfg
    python ..\plot_linear.py $cfg
    julia .\third_order_xcf.jl $cfg
    ls ${root}_r[nr].bin | %{julia .\spec_ft.jl $cfg $_ $_.Name.Replace("_r", "_s")}
    julia .\spec_abs.jl $cfg $root ${root}_sa.bin
    python ..\plot_2d.py $cfg
}
