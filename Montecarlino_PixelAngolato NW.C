#include "includefile.h"
Double_t GaussianaN(Double_t *x, Double_t *par){
  Double_t Norm    = par[0];
  Double_t mu    = par[1];
  Double_t sigma = par[2];
  return Norm/sqrt(6.28)/sigma*TMath::Exp(-pow((x[0]-mu),2)/(2.*pow(sigma,2)));
}
//=============================
Double_t GaussianaN(Double_t *x, Double_t *par);
Double_t BifurGauss(Double_t *x, Double_t *par);
//=============================

void Montecarlino_PixelAngolato(  Int_t Nwave=1000000){
  //gStyle->SetOptStat("");
  gStyle->SetOptStat("");
  gStyle->SetOptFit(1112);
  //gStyle->SetStatX(0.4);
  //gStyle->SetStatY(0.9);
  TRandom3 rndm;

  Float_t pix_depth =150;
  Float_t pix_size = 55.-10.;
  Float_t z_in = 50;
  Float_t x_in, x_out;
  Float_t z_out = z_in+pix_depth;
  Float_t _x_in, _x_out, _z_in,_z_out;
  Float_t max = 1.2*sqrt(pix_depth*pix_depth+pix_size*pix_size);
  Float_t min(0.);
  Int_t nbin=200;
  Float_t rot_x[6] = {0,3,5,10,15,20};
  TH1F *h_length[6], *h_ampl[6], *hX_in[6];
  TH1F *h_theta = new TH1F("h_theta","theta",50,-0.1,0.1);
  TH2F *h2D_in[6], *h2D_out[6];
  Float_t eff[6];
  Int_t icol[6]={1,2,4,6,7,9};
  TString hname("");
  TString htitle("");
  TLegend *legend = new TLegend(0.6-0.4,0.65,0.9-0.4,0.85);
  legend->SetLineColor(0);
  legend->SetFillStyle(0);

  // Landau parameters
  Float_t LandauRatio = 2.;           // FWHM_Landau/Peak_Landau
  Float_t Peak_at0degrees = 100.;     // fitted value for signal amplitude [mV]

  for(int i=0; i<6; i++){
    hname.Form("length%d",i);
    htitle.Form("angle%1.1f ",rot_x[i]);
    rot_x[i] = -rot_x[i]*3.14/180.;
    h_length[i] = new TH1F(hname, htitle,nbin,min,max);
    h_length[i]->SetLineColor(icol[i]);
    h_length[i]->SetLineWidth(2);
    h_length[i]->GetXaxis()->SetTitle("length [#mum]");
    hname.Form("h_xIN%d",i);
    htitle.Form("h_xIN%1.1f ",rot_x[i]);
    hX_in[i] = new TH1F(hname, htitle,nbin,-6.*pix_size,6*pix_size);
    hname.Form("IN%d",i);
    htitle.Form("IN%1.1f ",rot_x[i]);
    h2D_in[i] = new TH2F(hname, htitle,nbin,-pix_size,pix_size,nbin,z_in*0.8, z_out*1.2);
    hname.Form("OUT%d",i);
    htitle.Form("OUT%1.1f ",rot_x[i]);
    h2D_out[i] = new TH2F(hname, htitle,nbin,-pix_size,pix_size,nbin,z_in*0.8, z_out*1.2);
    for(int iw=0; iw<Nwave; iw++){
      Float_t weight=1;
      Float_t x = rndm.Uniform(-4*pix_size,4*pix_size);
      Float_t theta_x = rot_x[i]+ rndm.Gaus(0.,0.04);
      h_theta->Fill(theta_x);
      x_in  = x + tan(theta_x)*z_in;
      x_out = x + tan(theta_x)*z_out;
      hX_in[i]->Fill(x_in);
      /*
      cout<<
    " x_in="<<x_in<<" x_out="<<x_out<<
    " z_in="<<z_in<<" z_out="<<z_out<<endl;
      */
      if(fabs(x_in)<=pix_size/2.) {
    _z_in = z_in;
    _x_in = x_in;
      } else {
    if(x_in>pix_size/2. ){
      _x_in = pix_size/2.;
      _z_in = (pix_size/2. -x) /tan(theta_x);
      //cout<<"PIU:x="<<x<<"  x_in="<<_x_in<<" z_in="<<_z_in<<" calcolo="<<x+tan(theta_x)*_z_in<<" "<<theta_x<<endl;
    } else {
      _x_in = -pix_size/2.;
      _z_in = (-pix_size/2. - x) /tan(theta_x);
      //cout<<"MENO:x="<<x<<"  x_in="<<_x_in<<" z_in="<<_z_in<<" calcolo="<<x+tan(theta_x)*_z_in<<" "<<theta_x<<endl;
    }
    if(_z_in<z_in || _z_in>z_out ) continue;
      }
      if(fabs(x_out)<=pix_size/2.) {
    _z_out = z_out;
    _x_out = x_out;
      } else {
    if(x_out>pix_size/2. ) {
      _x_out = pix_size/2.;
      _z_out = (pix_size/2. -x) /tan(theta_x);
    } else {
      _x_out = -pix_size/2.;
      _z_out = (-pix_size/2. -x) /tan(theta_x);
    }
    if(_z_in<z_in || _z_in>z_out ) continue;
      }
      /*
      cout<<
    " x_in="<<x_in<<" x_out="<<x_out<<
    " z_in="<<z_in<<" z_out="<<z_out<<endl;
      */
      /*
      if(x_in>pix_size/2.&&theta_x<0) {
    _z_in = (pix_size/2.-x)/(x_in-x)*z_in;
    _x_in = pix_size/2.;
    if(_z_in<z_in || _z_in>z_out ) continue;
      }else if(x_in<-pix_size/2.&&theta_x>0) {
    _z_in = (-pix_size/2.-x)/(x_in-x)*z_in;
    _x_in = -pix_size/2.;
    if(_z_in<z_in || _z_in>z_out ) continue;
      }else if(fabs(x_in)<pix_size/2.) {
    _z_in = z_in;
    _x_in = x_in;
      } else {
    continue;
      }

      if(x_out>pix_size/2. &&theta_x>0 ) {
    _z_out = (pix_size/2.-x)/(x_out-x)*z_out;
    _x_out = pix_size/2.;
    if(_z_out<z_in || _z_out>z_out ) continue;

    //weight = 2.;
      }else if(x_out<-pix_size/2.&&theta_x<0) {
    _z_out = (-pix_size/2.-x)/(x_out-x)*z_out;
    _x_out = -pix_size/2.;
    if(_z_out<z_in || _z_out>z_out ) continue;
    //    weight = 2.;
      }else if(fabs(x_out)<pix_size/2.) {
    _z_out = z_out;
    _x_out = x_out;
      } else {
    continue;
      }
      */
      //      if(fabs(x_out)<pix_size/2.) continue;

      Float_t length = sqrt((_x_out-_x_in)*(_x_out-_x_in)+(_z_out-_z_in)*(_z_out-_z_in));

      Float_t scaled_peak = Peak_at0degrees*length/pix_depth;
      Float_t amplitude = rndm.Landau(scaled_peak,LandauRatio*scaled_peak);
      h_ampl[i]->Fill(amplitude);
      /*
      cout<<
    " _x_in="<<_x_in<<" _x_out="<<_x_out<<
    " _z_in="<<_z_in<<" _z_out="<<_z_out<<
    " length="<<length<<" | "<<(_x_out-_x_in)*(_x_out-_x_in)+(_z_out-_z_in)*(_z_out-_z_in)<<endl;
      */
      h2D_in[i]->Fill(_x_in,_z_in);
      h2D_out[i]->Fill(_x_out,_z_out);
      h_length[i]->Fill(length);
    }
}
  TCanvas *c_dist = new TCanvas("c","c",1000,500);
  c_dist->cd()->SetLogy();
  TString opt("");
  for(int i=0; i<6; i++){
    h_length[i]->Draw(opt);
    opt="same";
    eff[i] = h_length[i]->GetEntries()/float(Nwave);
    cout<<" Geometrical efficiency "<<htitle.Data()<<" "<<eff[i]<<endl;
    htitle.Form("#epsilon = %1.3f for angle %1.1f ",eff[i]/eff[0],rot_x[i]*180/3.14);
    legend->AddEntry(h_length[i],htitle,"l");
  }
  h_length[0]->Draw(opt);
  legend->Draw();

  TCanvas *c_ampl = new TCanvas("c_ampl","c_ampl",800,600);
    c_ampl->cd();
  for(int i=0; i<6; i++){
      h_ampl[i]->Draw(opt);
  }
//  h_ampl[0]->Draw(opt);

  //c2D->cd(1);  hX_in[0]->Draw();
  //c2D->cd(2);  h_theta->Draw();
  TCanvas *c2D =new TCanvas("c2D","c2D",900,900); c2D->Divide(6,2);
  for(int i=0; i<6; i++){
    c2D->cd(i+1);  h2D_in[i]->Draw("colz");
    c2D->cd(i+7);  h2D_out[i]->Draw("colz");
  }
}
