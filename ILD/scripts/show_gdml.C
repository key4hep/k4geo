
void show_gdml(const char* FILEN) {


TGeoManager::Import( FILEN );

gGeoManager->GetTopVolume()->Draw("ogl");

}
