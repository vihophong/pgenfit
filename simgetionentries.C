void simgetionentries()
{
TFile* f=TFile::Open("simoutput.root");
TTree* treeion = (TTree*) f->Get("ion");
TTree* treebeta = (TTree*) f->Get("beta");
std::ofstream ofs("simoutput_fit.root.txt",std::ios::app);
Long64_t nbeta = treebeta->Draw("","beta.id==0");
ofs<<treeion->GetEntries()<<"\t"<<nbeta<<endl;
ofs.close();
f->Close();
}
