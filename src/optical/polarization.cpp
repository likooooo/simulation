  // set PolarizationAngleMap                                                                      │
│   822             if (fPolarizationAngleFunc_) {                                                                   │
│   823                 MEM_LEAK_FLAG( oPolarizationAngleMap_ );                                                     │
│   824                 oPolarizationAngleMap_ = new CPageT<OLreal>;                                                 │
│   825                 oPolarizationAngleMap_->SetPageInit(spatialCoordsPixel_, fourierCoordsPixel_, true);         │
│   826                 for (OLint iy = 0; iy < oPolarizationAngleMap_->SizeY(); iy++) {                             │
│   827                     OLreal sigmaY = iy * oPolarizationAngleMap_->rStepY() + oPolarizationAngleMap_->rStartY()│
│   828                     for (OLint ix = 0; ix < oPolarizationAngleMap_->SizeX(); ix++) {                         │
│   829                         OLreal sigmaX = ix * oPolarizationAngleMap_->rStepX() + oPolarizationAngleMap_->rStar│
│  >830                         OLreal phi = fPolarizationAngleFunc_(sigmaX, sigmaY);                                │
│   831                         oPolarizationAngleMap_->SetOneWord(ix, iy, phi);                                     │
│   832                     }                                                                                        │
│   833                 }                                                                                            │
│   834             }         

sourcePolarization = "XY"
linearPolarization = 10
def SourcePolAngleFunc(fx, fy):
    if sourcePolarization == "XY":
        numFault = 0.00001 # Avoid the influence of numerical error
        if( ((fy-fx<-numFault) and (fy+fx>numFault)) or ((fy-fx>numFault) and (fy+fx<-numFault)) ):
            pol_angle = 90
        elif( ((fy-fx>numFault) and (fy+fx>numFault)) or ((fy-fx<-numFault) and (fy+fx<-numFault)) ):
            pol_angle = 0
        else : pol_angle = -361  # Special angle(FLAG) used to represent natural light
    if sourcePolarization == "X":
        pol_angle = 0
    if sourcePolarization == "Y":
        pol_angle = 90
    if sourcePolarization == "TE":
        numFault = 0.00001 # Avoid the influence of numerical error
        if ((abs(fy)<numFault) and (abs(fx)<numFault)):
            pol_angle = 90
        else:
            pol_angle = atan2(fy,fx)/pi*180+90        
    if sourcePolarization == "TM":
        numFault = 0.00001 # Avoid the influence of numerical error
        if ((abs(fy)<numFault) and (abs(fx)<numFault)):
            pol_angle = 0
        else:
            pol_angle = atan2(fy,fx)/pi*180         
    if sourcePolarization == "LP":
        pol_angle = linearPolarization
    if sourcePolarization == "UNPOLARIZED":
        pol_angle = -361
    return pol_angle # in units of degree
