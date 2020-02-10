$cfg=$args
if (!$cfg) {
    $cfg=$(ls *.yaml)
}
foreach ($c in $cfg) {
    $root=$(Select-String -Path $c -Pattern '^rootname:\W*(.*)$').Matches.Groups[1].value
    julia .\linear_xcf.jl $c
    python ..\plot_linear.py $c
    julia .\third_order_xcf.jl $c
    ls ${root}_r[nr].bin | %{julia .\spec_ft.jl $c $_ $_.Name.Replace("_r", "_s")}
    julia .\spec_abs.jl $c $root ${root}_sa.bin
    python ..\plot_2d.py $c
}
