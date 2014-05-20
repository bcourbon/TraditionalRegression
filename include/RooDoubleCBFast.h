#ifndef ROODOUBLECBFAST
#define ROODOUBLECBFAST

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"

class RooDoubleCBFast : public RooAbsPdf {
public:
  RooDoubleCBFast();
  RooDoubleCBFast(const char *name, const char *title,
              RooAbsReal& _x,
              RooAbsReal& _mean,
              RooAbsReal& _width,
              RooAbsReal& _alpha1,
              RooAbsReal& _n1,
              RooAbsReal& _alpha2,
              RooAbsReal& _n2
           );
  RooDoubleCBFast(const RooDoubleCBFast& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooDoubleCBFast(*this,newname); }
  inline virtual ~RooDoubleCBFast() { }

protected:

  RooRealProxy x ;
  RooRealProxy mean;
  RooRealProxy width;
  RooRealProxy alpha1;
  RooRealProxy n1;
  RooRealProxy alpha2;
  RooRealProxy n2;
  
  Double_t evaluate() const ;

private:

  ClassDef(RooDoubleCBFast,1)
};
#endif
