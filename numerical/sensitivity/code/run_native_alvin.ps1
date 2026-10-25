# Compile the two unmodified numerical methods from the pinned upstream C# source.
# This exercises the external pressure kernel in its original float arithmetic;
# Unity rendering and the Dash Hb module are outside this qualification.
$ErrorActionPreference = 'Stop'
$phase = Split-Path -Parent $PSScriptRoot
$source = Get-Content -LiteralPath (Join-Path $phase 'external_alvin/PartialPressure.cs') -Raw
function Extract-Method([string]$methodName) {
    $start = $source.IndexOf('private void ' + $methodName + '()')
    if ($start -lt 0) { throw "Missing upstream method $methodName" }
    $open = $source.IndexOf('{', $start)
    $depth = 1
    $end = $open + 1
    while ($depth -gt 0) {
        if ($source[$end] -eq '{') { $depth++ }
        if ($source[$end] -eq '}') { $depth-- }
        $end++
    }
    return $source.Substring($start, $end-$start)
}
$methods = (Extract-Method 'CharacterizeSections') + "`n" + (Extract-Method 'CalculatePressures')
$wrapper = @'
using System;
public class AlvinPressureKernel {
private int numberSections=50;
private const int MeanCapLength=500;
private const float volumeCapillaryNet=808000;
private const float RQ=.82f;
private float bloodVelocity,surfaceArea,tissueBarrier,pO2_alv,pO2_blood,pCO2_alv,pCO2_blood,atmosphericPressure;
private float sectionLength,timePeriod,sectionVolume,sectionArea;
private float[] partialPressures,o2GainPerSection,co2partialPressures,co2GainPerSection;
public float[] Run(float pv,float q,float beta,int sections) {
 numberSections=sections; pO2_alv=100; pO2_blood=pv; pCO2_alv=40;pCO2_blood=46;
 atmosphericPressure=760;tissueBarrier=1;
 float transit=.06f*86f/q;
 bloodVelocity=MeanCapLength/transit;
 surfaceArea=beta*volumeCapillaryNet*tissueBarrier/(760f*.055f);
 CharacterizeSections();CalculatePressures();return partialPressures;
}
'@
Add-Type -TypeDefinition ($wrapper + $methods + "`n}")
$summary = Get-Content -LiteralPath (Join-Path $phase 'results/run_summary.json') -Raw | ConvertFrom-Json
$beta=[float]$summary.benchmark_beta
$q0=[float]3.9464411557434818
$states=@(@(40,$q0),@(30,$q0),@(50,$q0),@(40,3),@(40,5))
$outputRows=foreach ($state in $states) {
 $kernel=New-Object AlvinPressureKernel
 $values=$kernel.Run([float]$state[0],[float]$state[1],$beta,50)
 for ($idx=0;$idx -lt $values.Length;$idx++) {
  [pscustomobject]@{Pv=$state[0];Q=[double]$state[1];node=$idx;native_PO2=[double]$values[$idx]}
 }
}
$outputRows | Export-Csv -LiteralPath (Join-Path $phase 'results/external_alvin_native.csv') -NoTypeInformation -Encoding utf8
Write-Output 'Compiled upstream CharacterizeSections and CalculatePressures; wrote 250 native C# node outputs.'
