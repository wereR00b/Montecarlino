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

void Montecarlino_PixelAngolato(Int_t Nwave=10000000){
    //gStyle->SetOptStat("");
  gStyle->SetOptStat("");
  gStyle->SetOptFit(1112);
    //gStyle->SetStatX(0.4);
    //gStyle->SetStatY(0.9);
  TRandom3 rndm;
  /*
   for(int i=0; i<6; i++){
   }
   */

  Float_t pix_depth =150.;  // !!! µm
  Float_t pix_size = 55.-10.;  // !!! µm
  Float_t trench_w = 5.;  // !!! µm
  Float_t z_in = 50.;  // !!! µm
  Float_t x_in, x_out;
  Float_t z_out = z_in+pix_depth;
  Float_t _x_in, _x_out, _z_in,_z_out;
  Float_t max = 1.2*sqrt(pix_depth*pix_depth+pix_size*pix_size);
  Float_t min(0.);
  Int_t nbin=1000;
  Float_t rotD_x[] = {0,3,10,20};  // rotation angles in degrees
  Int_t anglesN = sizeof(rotD_x)/sizeof(rotD_x[0]);
  Float_t rot_x[anglesN];
  TH1F *h_length[anglesN], *h_ampl[anglesN], *hX_in[anglesN];
  TH1F *h_theta = new TH1F("h_theta","theta",50,-0.1,0.1);
  TH2F *h2D_in[anglesN], *h2D_out[anglesN];
  Float_t eff[anglesN];
  Int_t icol[]={1,2,4,6,7,9};
  TString hname("");
  TString htitle("");
    // Landau parameters
  Float_t LandauRatio = 5.;          // LandauPeak/LandauFWHM
  Float_t LandauPeak_meas = 70.;    // fitted value for signal amplitude [mV]
  Float_t elNoise = 1.;             // sigma (RMS) value of electronic noise (Gaussian)

  for(int i=0; i<anglesN; i++){
    hname.Form("length%d",i);
    htitle.Form("Track length distribution for different tilt angles");
    rot_x[i] = rotD_x[i]*3.14/180.;
    h_length[i] = new TH1F(hname, htitle,nbin,min,max);
    h_length[i]->SetLineColor(icol[i]);
    h_length[i]->SetLineWidth(2);
    h_length[i]->GetXaxis()->SetTitle("length [#mum]");

    hname.Form("Amplitude%d",i);
    htitle.Form("Amplitude distribution for different tilt angles");
    //htitle.Form("Amplitude distribution for different tilt angles %1.0f",rot_x[i]);
    h_ampl[i] = new TH1F(hname, htitle,nbin,0,400);
    h_ampl[i]->SetLineColor(icol[i]);
    h_ampl[i]->SetLineWidth(2);
    h_ampl[i]->GetXaxis()->SetTitle("amplitude [mV]");
    hname.Form("h_xIN%d",i);
    htitle.Form("h_xIN%1.0f deg",rotD_x[i]);
    hX_in[i] = new TH1F(hname, htitle,nbin,-6.*pix_size,6*pix_size);
    hname.Form("IN%d",i);
    htitle.Form("IN %1.0f deg",rotD_x[i]);
    h2D_in[i] = new TH2F(hname, htitle,nbin,-pix_size,pix_size,nbin,z_in*0.8, z_out*1.2);
    hname.Form("OUT%d",i);
    htitle.Form("OUT %1.0f deg",rotD_x[i]);
    h2D_out[i] = new TH2F(hname, htitle,nbin,-pix_size,pix_size,nbin,z_in*0.8, z_out*1.2);
    for(int iw=0; iw<Nwave; iw++){
      Float_t weight=1;
      Float_t theta_x = rot_x[i]+ rndm.Gaus(0.,0.01); // sigma of the beam ~0.0001rad
      h_theta->Fill(theta_x);
      Float_t x = rndm.Uniform(-4*pix_size+z_in*tan(theta_x),4*pix_size+z_in*tan(theta_x)); // x coordinate of the vertex (z=0)
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
      Float_t length = sqrt((_x_out-_x_in)*(_x_out-_x_in)+(_z_out-_z_in)*(_z_out-_z_in));
      if (fabs(_x_in)<0.5*trench_w && fabs(_x_out)<0.5*trench_w) {  // whole track is inside trench
        continue;
      } else if (fabs(_x_in)<0.5*trench_w){
     // Float_t DeltaW = 0.5*trench_w+ theta_x/fabs(theta_x)*_x_in;   //-> trench width intercepted by track
        length = length - (0.5*trench_w+ theta_x/fabs(theta_x)*_x_in)/fabs(sin(theta_x));
      } else if (fabs(_x_out)<0.5*trench_w){
     // Float_t DeltaW = 0.5*trench_w- theta_x/fabs(theta_x)*_x_out;   //-> trench width intercepted by track
        length = length - (0.5*trench_w- theta_x/fabs(theta_x)*_x_out)/fabs(sin(theta_x));
      } else if (_x_in*_x_out<0.){
        length = length - trench_w/fabs(sin(theta_x));
      }

      if (length<0.) continue; // select track longer than a threshold

      Float_t scaled_peak = LandauPeak_meas*length/pix_depth;
      Float_t amplitude = rndm.Landau(scaled_peak,scaled_peak/LandauRatio)+rndm.Gaus(0.,elNoise);
        //        cout<<"length= "<<length<<" scaled_peak= "<<scaled_peak<<" amplitude= "<<amplitude<<endl;
      h_ampl[i]->Fill(amplitude);
      /*
       cout<< " _x_in="<<_x_in<<" _x_out="<<_x_out<<" _z_in="<<_z_in<<" _z_out="<<_z_out<<
       " length="<<length<<" | "<<(_x_out-_x_in)*(_x_out-_x_in)+(_z_out-_z_in)*(_z_out-_z_in)<<endl;
       */
      h2D_in[i]->Fill(_x_in,_z_in);
      h2D_out[i]->Fill(_x_out,_z_out);
      h_length[i]->Fill(length);
    }
  }
  TCanvas *c = new TCanvas("c","c",900,450); c->Divide(2,1);
  TLegend *legend = new TLegend(0.15,0.4,0.5,0.2);
  legend->SetLineColor(0);
  legend->SetFillStyle(0);

  c->cd(1)->SetLogy();
  TString opt("");
  for(int i=0; i<anglesN; i++){
    h_length[i]->Draw(opt);
    opt="same";
    eff[i] = h_length[i]->GetEntries()/float(Nwave);
    //cout<<" Geometrical efficiency "<<htitle.Data()<<" "<<eff[i]<<endl;
    htitle.Form("#epsilon = %1.2f for angle %1.0f deg ",eff[i]/eff[0],rot_x[i]*180/3.14);
    legend->AddEntry(h_length[i],htitle,"l");
  }
  h_length[0]->Draw(opt);
  legend->Draw();
  c->cd(2);

  TLegend *legend2 = new TLegend(0.6-0.1,0.65,0.9-0.1,0.85);
  legend2->SetLineColor(0);
  legend2->SetFillStyle(0);
  for(int i=anglesN-1; i>-1; i--){
    opt="same";
    h_ampl[i]->Draw(opt);
    htitle.Form("angle %1.0f deg ",rot_x[i]*180/3.14);
    legend2->AddEntry(h_length[i],htitle,"l");
  }
  legend2->Draw();

    //c2D->cd(1);  hX_in[0]->Draw();
    //c2D->cd(2);  h_theta->Draw();
  TCanvas *c2D =new TCanvas("c2D","c2D",anglesN*300,800); c2D->Divide(anglesN,2);
  for(int i=0; i<anglesN; i++){
    c2D->cd(i+1);  h2D_in[i]->Draw("colz");
    c2D->cd(i+anglesN+1);  h2D_out[i]->Draw("colz");
  }
}
