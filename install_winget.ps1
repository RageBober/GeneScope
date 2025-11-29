# Download and install winget
$tempPath = Join-Path $env:TEMP "winget.msixbundle"
$vcLibsPath = Join-Path $env:TEMP "Microsoft.VCLibs.x64.14.00.Desktop.appx"

Write-Host "Downloading winget..."
curl.exe -L "https://aka.ms/getwinget" -o $tempPath

Write-Host "Downloading VC++ Runtime..."
curl.exe -L "https://aka.ms/Microsoft.VCLibs.x64.14.00.Desktop.appx" -o $vcLibsPath

Write-Host "Installing VC++ Runtime..."
Add-AppxPackage -Path $vcLibsPath

Write-Host "Installing winget..."
Add-AppxPackage -Path $tempPath

Write-Host "Done! Cleaning up..."
Remove-Item $tempPath -ErrorAction SilentlyContinue
Remove-Item $vcLibsPath -ErrorAction SilentlyContinue

Write-Host "Installation complete!"
