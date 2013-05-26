void Dopen()
{
	TCanvas *c = new TCanvas("c","omega",0,0,700,600);
	TGraph *dt = new TGraph("/home/alexandra/CP_Hw/Chap3_2/FFT==>>4-order_RK_Method.dat");
//	gStyle->SetPalette(1);
	dt->Draw("AP");
//	return c;
}

