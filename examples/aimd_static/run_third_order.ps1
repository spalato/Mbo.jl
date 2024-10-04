$cfg=$args
if (!$cfg) {
    $cfg=$(ls *.yaml)
}
foreach ($c in $cfg) {
    
    $root=$(Select-String -Path $c -Pattern '^rootname:\W*(.*)$').Matches.Groups[1].value
    julia --project .\third_order_xc.jl $c
    ls ${root}_r[nr].bin | %{julia --project .\spec_ft.jl $c $_ $_.Name.Replace("_r", "_s")}
    ls ${root}_s[nr].bin | %{python .\plot_spec.py $c $_}
    julia --project .\spec_abs.jl $c $root ${root}_sa.bin
    #rm *_[sr][rn].bin
    python .\plot_spec.py $c ${root}_sa.bin
    #python .\beatmap.py $cfg ${root}_sa.bin
}