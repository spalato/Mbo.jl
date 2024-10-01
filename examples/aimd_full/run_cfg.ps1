$cfg=$args
if (!$cfg) {
    $cfg=$(ls *.yaml)
}
foreach ($c in $cfg) {
    $root=$(Select-String -Path $c -Pattern '^rootname:\W*(.*)$').Matches.Groups[1].value
    julia --project .\linear_xcf.jl $c
    python ..\plot_linear.py $c
    julia --project .\third_order_xcf.jl $c
    ls ${root}_r[nr].bin | %{julia --project .\spec_ft.jl $c $_ $_.Name.Replace("_r", "_s")}
    julia --project .\spec_abs.jl $c $root ${root}_sa.bin
    python ..\plot_2d.py $c
}
