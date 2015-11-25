/*
-------------------------------------------------
-------     Universidad de los Andes      -------
-------      Departamento de Física       -------
-------        Joven Investigador         -------
-------  Andrés Felipe García Albarracín  -------
-------    Juan Carlos Sanabria Arenas    -------
-------------------------------------------------

This algorithm studies the kinematic properties
of the ISR jets. It reads the results of the
matching algorithm

To execute, type:

./ISR_jet_analysis config_file.txt

where config_file.txt is the mandatory configuration
file
*/


#include "ROOTFunctions.h"
#include "graphs_Funcs.h"
#include "functions.h"
#include "Rtypes.h"
#include "DelphesFunctions.h"
#include "mt2_bisect.h"
#include "mt2w_bisect.h"
#include "mt2bl_bisect.h"



// Global Variables
const Double_t PI = TMath::Pi();

float mtop=173.34;
float mw=80.40;
float sigma_12= 10.5;
float sigma_123 = 25.7;
double mn    = 0.;
double pi = atan(1.0)*4.0;



int main(int argc, char **argv){
  std::cout.precision(4);
  // Counting time
  Double_t initialTime = clock();
  
  // Folder variables
  // Masa del stop de 200 GeV y del LSP de 20 GeV
         
  string head_folder = "/home/rrodriguez/Simulación/MG_pythia8_delphes_parallel/Runs/";
  string current_folder = "Susy_stops1_13Tev_con_1_ISR_stop_200_N1_20_001/";
  string head_folder_binary = "/home/rrodriguez/Simulación/Results_Improved_Codes_Susy_stops_200_N1_20/matching_Results/Stops_200_N1_20_matchs_WI_Matching/";
  string matching_name = "ISR_jets_Susy_stops_200_N1_20_WI_001.bn";
  string head_folder_results =  "/home/rrodriguez/Simulación/Results_Improved_Codes_Susy_stops_200_N1_20/resultsTagging/Stops_200_N1_20_histos_WI_Matching/";
  
  //Masa del stop de 200 GeV y del LSP de 35 GeV
  // Folder variables
  /*         
  string head_folder = "/home/rrodriguez/Simulación/MG_pythia8_delphes_parallel/Runs/";
  string current_folder = "Susy_stops1_13Tev_con_1_ISR_stop_200_N1_35_001/";
  string head_folder_binary = "/home/rrodriguez/Simulación/Results_Improved_Codes_Susy_stops_200_N1_35/matching_Results/Stops_200_N1_35_matchs_WI_Matching/";
  string matching_name = "ISR_jets_Susy_stops_200_N1_35_WI_001.bn";
  string head_folder_results =  "/home/rrodriguez/Simulación/Results_Improved_Codes_Susy_stops_200_N1_35/resultsTagging/Stops_200_N1_35_histos_WI_Matching/";
  */
 //Masa del stop de 300 GeV y del LSP de 135 GeV
  // Folder variables
  /*              
  string head_folder = "/home/rrodriguez/Simulación/MG_pythia8_delphes_parallel/Runs/";
  string current_folder = "Susy_stops1_13Tev_con_1_ISR_stop_300_N1_135_001/";
  string head_folder_binary = "/home/rrodriguez/Simulación/Results_Improved_Codes_Susy_stops_300_N1_135/matching_Results/Stops_300_N1_135_matchs_WI_Matching/";
  string matching_name = "ISR_jets_Susy_stops_300_N1_135_WI_001.bn";
  string head_folder_results =  "/home/rrodriguez/Simulación/Results_Improved_Codes_Susy_stops_300_N1_135/resultsTagging/Stops_300_N1_135_histos_WI_Matching/";
  */
 //Masa del stop de 300 GeV y del LSP de 120 GeV
  // Folder variables
  /*             
  string head_folder = "/home/rrodriguez/Simulación/MG_pythia8_delphes_parallel/Runs/";
  string current_folder = "Susy_stops1_13Tev_con_1_ISR_stop_300_N1_120_001/";
  string head_folder_binary = "/home/rrodriguez/Simulación/Results_Improved_Codes_Susy_stops_300_N1_120/matching_Results/Stops_300_N1_120_matchs_WI_Matching/";
  string matching_name = "ISR_jets_Susy_stops_300_N1_120_WI_001.bn";
  string head_folder_results =  "/home/rrodriguez/Simulación/Results_Improved_Codes_Susy_stops_300_N1_120/resultsTagging/Stops_300_N1_120_histos_WI_Matching/";
  */
  //Masa del stop de 400 GeV y del LSP de 220 GeV
  // Folder variables
  /*           
  string head_folder = "/home/rrodriguez/Simulación/MG_pythia8_delphes_parallel/Runs/";
  string current_folder = "Susy_stops1_13Tev_con_1_ISR_stop_400_N1_220_001/";
  string head_folder_binary = "/home/rrodriguez/Simulación/Results_Improved_Codes_Susy_stops_400_N1_220/matching_Results/Stops_400_N1_220_matchs_WI_Matching/";
  string matching_name = "ISR_jets_Susy_stops_400_N1_220_WI_001.bn";
  string head_folder_results =  "/home/rrodriguez/Simulación/Results_Improved_Codes_Susy_stops_400_N1_220/resultsTagging/Stops_400_N1_220_histos_WI_Matching/";
  */
  //Masa del stop de 400 GeV y del LSP de 235 GeV
  // Folder variables
  /*            
  string head_folder = "/home/rrodriguez/Simulación/MG_pythia8_delphes_parallel/Runs/";
  string current_folder = "Susy_stops1_13Tev_con_1_ISR_stop_400_N1_235_001/";
  string head_folder_binary = "/home/rrodriguez/Simulación/Results_Improved_Codes_Susy_stops_400_N1_235/matching_Results/Stops_400_N1_235_matchs_WI_Matching/";
  string matching_name = "ISR_jets_Susy_stops_400_N1_235_WI_001.bn";
  string head_folder_results =  "/home/rrodriguez/Simulación/Results_Improved_Codes_Susy_stops_400_N1_235/resultsTagging/Stops_400_N1_235_histos_WI_Matching/";
  */
  //Masa del stop de 500 GeV y del LSP de 320 GeV
  // Folder variables
  /*          
  string head_folder = "/home/rrodriguez/Simulación/MG_pythia8_delphes_parallel/Runs/";
  string current_folder = "Susy_stops1_13Tev_con_1_ISR_stop_500_N1_320_001/";
  string head_folder_binary = "/home/rrodriguez/Simulación/Results_Improved_Codes_Susy_stops_500_N1_320/matching_Results/Stops_500_N1_320_matchs_WI_Matching/";
  string matching_name = "ISR_jets_Susy_stops_500_N1_320_WI_001.bn";
  string head_folder_results =  "/home/rrodriguez/Simulación/Results_Improved_Codes_Susy_stops_500_N1_320/resultsTagging/Stops_500_N1_320_histos_WI_Matching/";
  */
  //Masa del stop de 400 GeV y del LSP de 100 GeV
  // Folder variables
  /*         
  string head_folder = "/home/rrodriguez/Simulación/MG_pythia8_delphes_parallel/Runs/";
  string current_folder = "Susy_stops1_13Tev_con_1_ISR_001/";
  string head_folder_binary = "/home/rrodriguez/Simulación/Results_Improved_Codes/matching_Results/Stops1_matchs_WI_Matching/";
  string matching_name = "ISR_jets_Stops1_WI_001.bn";
  string head_folder_results =  "/home/rrodriguez/Simulación/Results_Improved_Codes/resultsTagging/Stops1_histos_WI_Matching/";
  */
  //Background Semileptónico
 // Folder variables
  /*       
  string head_folder = "/Disco1/Pheno/BackgroudSamples/ttbar_semileptonico/";
  string current_folder = "Top_semileptónico_13TeV_con_ISR_001/";
  string head_folder_binary = "/home/rrodriguez/Simulación/Results_Improved_Codes_Top_semileptónico/matching_Results/Top_semileptónico_matchs_WI_Matching/";
  string matching_name = "ISR_jets_Top_semileptónico_1_WI_001.bn";
  string head_folder_results =  "/home/rrodriguez/Simulación/Results_Improved_Codes_Top_semileptónico/resultsTagging/Top_semileptónico_histos_WI_Matching/";
  */

  //Background Dileptónico
  // Folder variables
  /*           
  string head_folder = "/Disco1/Pheno/BackgroudSamples/ttbar_dileptonico/";
  string current_folder = "Top_leptónico_13TeV_con_ISR_001/";
  string head_folder_binary = "/home/rrodriguez/Simulación/Results_Improved_Codes_Top_leptónico/matching_Results/Top_leptónico_matchs_WI_Matching/";
  string matching_name = "ISR_jets_Top_leptónico_1_WI_001.bn";
  string head_folder_results =  "/home/rrodriguez/Simulación/Results_Improved_Codes_Top_leptónico/resultsTagging/Top_leptónico_histos_WI_Matching/";
  */
  // Checking input parameters
  string config_file_name = "Debug/config_file.txt";
  // Reading the file as first parameter
  if (argc>1){
    config_file_name = argv[1];
  }
  else{
    cout << "It is necessary to type a configuration file as parameter. Execute as ./ISR_jet_analysis config_file.txt" << endl;
    return 1;
  }
  cout << "Reading input parameters" << endl;
  cout << "\tUsing as parameters' file: " << config_file_name << endl;
  
  ifstream config_file (config_file_name);
  if (config_file.is_open()){
    cout << "\tReading file" << endl;
    string line;
    int number_line = 1;
		while (getline(config_file,line)){
		  // Skipping commented lines
		  if (line[0] == '!')
		    continue;
		  
		  // Finding the position of the equal sign
		  int pos_equal = -1;
		  pos_equal = line.find('=');
		  
		  if (pos_equal == -1){
		    cout << "\tLine " << number_line << " is incorrect" << endl;
		    continue;
		  }
		  
		  // Splitting the line according to the position of equal sign
		  string var_name = line.substr(0,pos_equal);
		  string var_value = line.substr(pos_equal+1);
		  
		  // Reading head folder
		  if(var_name.compare("head_folder") == 0){
		    head_folder = var_value;
		    cout << "\tVariable head folder set as: " << head_folder << endl;
		  }
		  // Reading current folder
		  else if (var_name.compare("current_folder") == 0){
		    current_folder = var_value;
		    cout << "\tVariable current folder set as: " << current_folder <<endl;
		  }
		  // Reading head folder binary
		  else if (var_name.compare("head_folder_binary") == 0){
		    head_folder_binary = var_value;
		    cout << "\tVariable head folder binary set as: " << head_folder_binary << endl;
		  }
		  // Reading matching name
		  else if (var_name.compare("matching_name") == 0){
		    matching_name = var_value;
		    cout << "\tVariable matching_name set as: " << matching_name << endl;
		  }
		  // Reading head folder results
		  else if (var_name.compare("head_folder_results") == 0){
		    head_folder_results = var_value;
		    cout << "\tVariable head folder results set as: " << head_folder_results << endl;
		  }
		  
		  number_line ++;
		}
  }
  else
    {
	    cout << "ERROR: File " << config_file_name << " does not exist. Terminating program" << endl;
	    return 0;
    }
  
  
	cout << "\n *** ISR jet analysis *** \n" << endl;
	
	/*
	 * Histograms
	 */
	// All jets
	
	TH1 *h_numberJet = new TH1F("Number Jets","Number Jets",11,-0.5,10.5);
	
	// Non Isr jets
	TH1 *h_jet_PT = new TH1F("Jet PT","Jet PT", 201,0.0,600.0);
	TH1 *h_jet_Eta = new TH1F("Jet Eta","Jet Eta", 171,-5.0,5.0);
	TH1 *h_jet_Phi = new TH1F("Jet Phi","Jet Phi", 375,-3.5,3.5);
	TH1 *h_jet_DPhi_MET = new TH1F("Jet - MET Delta_Phi","Jet - MET Delta_Phi",300,0.0,4.0);
	TH1 *h_jet_DPhi_MET_hpt = new TH1F("Jet - MET Delta_Phi_hpt","Jet - MET Delta_Phi_hpt",300,0.0,4.0);
	TH1 *h_jet_MT = new TH1F("Jet Transverse mass","Jet Transverse Mass",201,0.0,600.0);
	TH1 *h_jet_Delta_PT = new TH1F("Jet Delta-PT","Non ISR Delta-PT", 201,0.0,300.0);
	TH1 *h_jet_PT_HT = new TH1F("Jet PT-HT ratio","Jet PT-HT ratio",201,-0.0025,1.0025);
	TH1 *h_jet_PT_over_PT_others = new TH1F("Jet PT/PT_others","Jet PT/PT_others",401,-0.0025,2.0025);
	TH1 *h_jet_Eta_over_Eta_others = new TH1F("Jet Eta/Eta_others","Jet Eta/Eta_others",401,-0.0025,2.0025);
	TH1 *h_jet_DPhi_over_Phi_others = new TH1F("Jet Phi/Phi_others","Jet Phi/Phi_others",401,-0.0025,2.0025);
	TH1 *h_jet_Delta_Eta = new TH1F("Jet Delta-Eta","Jet Delta-Eta", 171,0.0,5.0);
	TH1 *h_jet_DPhi_MET_other = new TH1F("Jet - MET Delta_Phi other","Jet - MET Delta_Phi other",300,0.0,4.0);
	TH1 *h_jet_multiplicity = new TH1F("Jet - Multiplicity","Jet - Multiplicity",101,-0.5,100.5);
	TH1 *h_jet_DeltaR = new TH1F ("Jet - Delta_R","Jet - Delta_R",201,-0.0025,0.8025);
	TH1 *h_jet_Delta_PT_leading = new TH1F("Delta PT: leading - Jet","Delta PT: leading - Jet", 201,0.0,600.0);
	TH1 *h_jet_Delta_Eta_leading = new TH1F("Delta Eta: Jet - leading","Delta Eta: Jet - leading", 171,0.0,8.0);
	TH2 *h2_jet_PTEta=new TH2F("Non_ISR_Jet_PT_Eta","Non ISR Jet PT Vs. Eta",201,-1.25,501.25,201,-4.02,4.02);

	// ISR jets
	TH1 *h_ISR_PT = new TH1F("ISR PT","ISR PT", 201,0.0,600.0);
	TH1 *h_ISR_Eta = new TH1F("ISR Eta","ISR Eta", 171,-5.0,5.0);
	TH1 *h_ISR_Phi = new TH1F("ISR Phi","ISR Phi", 375,-3.5,3.5);
	TH1 *h_ISR_DPhi_MET = new TH1F("ISR - MET Delta_Phi","ISR - MET Delta_Phi",300,0.0,4.0);
	TH1 *h_ISR_DPhi_MET_hpt = new TH1F("ISR - MET Delta_Phi_hpt","ISR - MET Delta_Phi_hpt",300,0.0,4.0);
	//MT con corte en el PT del ISR
	TH1 *h_MT_Muon = new TH1F("ISR Transverse mass_Muon","ISR Transverse Mass muon",201,0.0,300.0);
	TH1 *h_MT_square_Muon = new TH1F("ISR Transverse mass_square_Muon","ISR Transverse Mass square muon",201,0.0,15000.0);
	TH1 *h_MT_Muon_corte_100_GeV = new TH1F("MT high_ISR_pt_corte_100_GeV_Muon","Transverse mass corte 100 GeV muon",200,0.01,250.0);
	TH1 *h_MT_Muon_corte_300_GeV = new TH1F("MT high_ISR_pt_corte_300_GeV_Muon","Transverse mass corte 300 GeV muon",200,0.01,250.0);
	TH1 *h_MT_Muon_corte_500_GeV = new TH1F("MT high_ISR_pt_corte_500_GeV_Muon","Transverse mass corte 500 GeV muon",200,0.01,250.0);
	TH1 *h_MT_Muon_square_corte_100_GeV = new TH1F("MT high_ISR_pt_square_corte_100_GeV_Muon","Transverse mass square corte 100 GeV muon",200,0.0,15000.0);
	TH1 *h_MT_Muon_square_corte_300_GeV = new TH1F("MT high_ISR_pt_square_corte_300_GeV_Muon","Transverse mass square corte 300 GeV muon",200,0.0,15000.0);
	TH1 *h_MT_Muon_square_corte_500_GeV = new TH1F("MT high_ISR_pt_square_corte_500_GeV_Muon","Transverse mass aquare corte 500 GeV muon",200,0.0,15000.0);
	TH1 *h_MT_Electron = new TH1F("ISR Transverse mass_Electron","ISR Transverse Mass electron",201,0.0,300.0);
	TH1 *h_MT_square_Electron = new TH1F("ISR Transverse mass_square_Electron","ISR Transverse Mass square electron",201,0.0,15000.0);
	TH1 *h_MT_Electron_corte_100_GeV = new TH1F("MT high_ISR_pt_corte_100_GeV_Electron","Transverse mass corte 100 GeV electron",200,0.01,250.0);
	TH1 *h_MT_Electron_corte_300_GeV = new TH1F("MT high_ISR_pt_corte_300_GeV_Electron","Transverse mass corte 300 GeV electron",200,0.01,250.0);
	TH1 *h_MT_Electron_corte_500_GeV = new TH1F("MT high_ISR_pt_corte_500_GeV_Electron","Transverse mass corte 500 GeV electron",200,0.01,250.0);
	TH1 *h_MT_Electron_square_corte_100_GeV = new TH1F("MT high_ISR_pt_square_corte_100_GeV_Electron","Transverse mass square corte 100 GeV electron",200,0.0,15000.0);
	TH1 *h_MT_Electron_square_corte_300_GeV = new TH1F("MT high_ISR_pt_square_corte_300_GeV_Electron","Transverse mass square corte 300 GeV electron",200,0.0,15000.0);
	TH1 *h_MT_Electron_square_corte_500_GeV = new TH1F("MT high_ISR_pt_square_corte_500_GeV_Electron","Transverse mass square corte 500 GeV electron",200,0.0,15000.0);
	//Mas histogramas
	TH1 *h_ISR_Delta_PT = new TH1F("ISR Delta-PT","ISR Delta-PT", 201,0.0,300.0);
	TH1 *h_ISR_PT_HT = new TH1F("ISR PT-HT ratio","ISR PT-HT ratio",201,-0.0025,1.0025);
	TH1 *h_ISR_PT_over_PT_others = new TH1F("ISR PT/PT_others","ISR PT/PT_others",401,-0.0025,2.0025);
	TH1 *h_ISR_Eta_over_Eta_others = new TH1F("ISR Eta/Eta_others","ISR Eta/Eta_others",401,-0.0025,2.0025);
	TH1 *h_ISR_DPhi_over_Phi_others = new TH1F("ISR Phi/Phi_others","ISR Phi/Phi_others",401,-0.0025,2.0025);
	TH1 *h_ISR_Delta_Eta = new TH1F("ISR Delta-Eta","ISR Delta-Eta", 171,0.0,5.0);
	TH1 *h_ISR_DPhi_MET_other = new TH1F("ISR - MET Delta_Phi other","ISR - MET Delta_Phi other",300,0.0,4.0);
	TH1 *h_ISR_multiplicity = new TH1F("ISR - Multiplicity","ISR - Multiplicity",101,-0.5,100.5);
	TH1 *h_ISR_DeltaR = new TH1F ("ISR - Delta_R","ISR - Delta_R",201,-0.0025,0.8025);
	TH1 *h_ISR_Delta_PT_leading = new TH1F("Delta PT: leading - ISR","Delta PT: leading - ISR", 201,0.0,600.0);
	TH1 *h_ISR_Delta_Eta_leading = new TH1F("Delta Eta: ISR - leading","Delta Eta: ISR - leading", 171,0.0,8.0);
	TH1 *h_ISR_MET_HT = new TH1F("ISR_MET_HT","isr_met_ht", 200,0.0,100);
	TH1 *h_ISR_MET_HT1 = new TH1F("ISR_MET_HT1","isr_met_ht1", 200,0.0,60);
	TH1 *h_ISR_MET_HT2 = new TH1F("ISR_MET_HT2","isr_met_ht1", 200,0.0,60);
	TH1 *h_ISR_MET_HT3 = new TH1F("ISR_MET_HT3","isr_met_ht1", 200,0.0,60);
	TH2 *h2_ISR_PTEta=new TH2F("ISR_Jet_PT_Eta","ISR Jet PT Vs. Eta",201,-1.25,501.25,201,-4.02,4.02);

	// MET con corte en el PT del ISR
	TH1 *h_MET = new TH1F("Missing ET","Missing ET",200,0,600);
	TH1 *h_MET_corte_100_GeV = new TH1F("Missing ET high_ISR_pt_corte_100_GeV","Missing ET corte 100 GeV",200,0.0,800.0);
	TH1 *h_MET_corte_300_GeV = new TH1F("Missing ET high_ISR_pt_corte_300_GeV","Missing ET corte 300 GeV",200,0.0,800.0);
	TH1 *h_MET_corte_500_GeV = new TH1F("Missing ET high_ISR_pt_corte_500_GeV","Missing ET corte 500 GeV",200,0.0,800.0);

	TH2 *h2_dif_PTEta=new TH2F("FSR_ISR_Jet_PT_Eta_Difference","Difference between FSR and ISR Jet PT Vs. Eta distributions",201,-1.25,501.25,201,-4.02,4.02);
	TH2 *h2_dif_lead_PTEta=new TH2F("Lead_ISR_Jet_PT_Eta_Difference","Difference between Lead and ISR Jet PT Vs. Eta distributions",201,-1.25,501.25,201,-4.02,4.02);

	// Leading PT
	TH1 *h_leading_PT = new TH1F("Leading PT","Leading PT", 201,0.0,600.0);
	TH1 *h_leading_MT = new TH1F("Leading Transverse mass","Leading Transverse Mass",201,0.0,600.0);
	TH1 *h_leading_Eta = new TH1F("Leading Eta","Leading Eta", 171,-5.0,5.0);
	TH1 *h_leading_DPhi_MET = new TH1F("Leading - MET Delta_Phi","Leading - MET Delta_Phi",300,0.0,4.0);

	TH2 *h2_leading_PTEta=new TH2F("Leading_Jet_PT_Eta","Leading Jet PT Vs. Eta",201,-1.25,501.25,201,-4.02,4.02);

	// Other variables
	TH1 *h_HT = new TH1F("HT","HT",201,0.0,600.0);
	TH1 *h_HT_R1 = new TH1F("HT_R1","HT_R1",51,-0.01,1.01);
	TH1 *h_HT_R2 = new TH1F("HT_R2","HT_R2",51,-0.01,1.01);

	// B tagging
	TH1 *h_BTag = new TH1F("BTag","BTag",5,-0.5,4.5);
	TH1 *h_BTag_PT = new TH1F("BTag PT","BTag PT", 201,0.0,600.0);
	TH1 *h_BTag_Eta = new TH1F("BTag Eta","BTag Eta", 171,-5.0,5.0);
	TH1 *h_BTag_DPhi_MET = new TH1F("BTag - MET Delta_Phi","BTag - MET Delta_Phi",300,0.0,4.0);
	TH1 *h_BTags_per_Event = new TH1F("BTags per event","BTags per event",5,-0.5,4.5);

	// Further analysis
	TH1 *h_ISR_PT_comp = new TH1F("ISR PT for comparison","ISR PT for comparison with histo", 20,0.0,800.0);
	TH1 *h_ISR_Eta_comp = new TH1F("ISR Eta for comparison","ISR Eta for comparison with histo", 20,-4.2,4.2);
	TH1 *h_ISR_DPhi_MET_comp = new TH1F("ISR Phi for comparison","ISR Phi for comparison with histo", 20,0,PI);

	// To check the histograms' creation
	TH1 *hist_ISR_PT = new TH1F("ISR PT comp","ISR PT comp", 20,0.0,800.0);
	TH1 *hist_ISR_Abs_Eta = new TH1F("ISR Abs Eta comp","ISR Abs Eta comp", 20,0.0,5.2);
	TH1 *hist_ISR_DPhi_MET = new TH1F("ISR Delta Phi comp","ISR Delta Phi comp", 20,0.0,PI);
	TH1 *hist_ISR_PT_ratio = new TH1F("ISR PT/PT_others comp","ISR PT/PT_others comp",20,0.0,8.0);
	TH1 *hist_ISR_Delta_Eta = new TH1F("ISR Delta-Eta comp","ISR Delta-Eta comp", 20,0.0,7.0);
	TH1 *hist_ISR_DPhi_MET_other = new TH1F("ISR - MET Delta_Phi other comp","ISR - MET Delta_Phi other comp",20,0.0,PI);
	TH1 *hist_ISR_Delta_PT_leading = new TH1F("Delta PT: leading - ISR comp","Delta PT: leading - ISR comp", 20,0.0,500.0);
	TH1 *hist_ISR_Delta_Eta_leading = new TH1F("Delta Eta: ISR - leading comp","Delta Eta: ISR - leading comp", 20,0.0,6.5);
	TH1 *hist_jet_PT = new TH1F("Jet PT comp","Jet PT comp", 20,0.0,800.0);
	TH1 *hist_jet_Abs_Eta = new TH1F("Jet Abs Eta comp","Jet Abs Eta comp", 20,0.0,5.2);
	TH1 *hist_jet_DPhi_MET = new TH1F("Jet Delta Phi comp","Jet Delta Phi comp", 20,0.0,PI);
	TH1 *hist_jet_PT_ratio = new TH1F("Jet PT/PT_others comp","Jet PT/PT_others comp",20,0.0,7.0);
	TH1 *hist_jet_Delta_Eta = new TH1F("Jet Delta-Eta comp","Jet Delta-Eta comp", 20,0.0,8.0);
	TH1 *hist_jet_DPhi_MET_other = new TH1F("Jet - MET Delta_Phi other comp","Jet - MET Delta_Phi other comp",20,0.0,PI);
	TH1 *hist_jet_Delta_PT_leading = new TH1F("Delta PT: leading - Jet comp","Delta PT: leading - Jet comp", 20,0.0,500.0);
	TH1 *hist_jet_Delta_Eta_leading = new TH1F("Delta Eta: Jet - leading comp","Delta Eta: Jet - leading comp", 20,0.0,6.5);

	//Histogramas de masa stransversa
	TH1 *histMT2 = new TH1F("MT2","mt2",100, 100.0, 490.0);
	TH1 *histMT2_corte_100_GeV = new TH1F("MT2_100_GeV","mt2_100_gev",100, 100.0, 490.0);
	TH1 *histMT2_corte_300_GeV = new TH1F("MT2_300_GeV","mt2_300_gev",100, 100.0, 490.0);
	TH1 *histMT2_corte_500_GeV = new TH1F("MT2_500_GeV","mt2_500_gev",100, 100.0, 490.0);
        TH1 *histMTB_Square = new TH1F("MTB_Square","mtb_square",100, 100.0, 25000.0);
        TH1 *histMTB_Square_corte_100_GeV = new TH1F("MTB_Square_100_GeV","mtb_square_100_GeV",100, 100.0, 25000.0);
        TH1 *histMTB_Square_corte_300_GeV = new TH1F("MTB_Square_300_GeV","mtb_square_300_GeV",100, 100.0, 25000.0);
        TH1 *histMTB_Square_corte_500_GeV = new TH1F("MTB_Square_500_GeV","mtb_square_500_GeV",100, 100.0, 25000.0);
        TH1 *hist_MT_square_MTB_square = new TH1F("MTB_Square_MT_square","mtb_square_MTB",100, 100.0, 30000.0);
	TH1 *histMT2W = new TH1F("MT2W","mt2w",100, 100.0, 490.0);
	TH1 *histMT2W_corte_100_GeV = new TH1F("MT2W_100_GeV","mt2w_100_gev",100, 100.0, 490.0);
	TH1 *histMT2W_corte_300_GeV = new TH1F("MT2W_300_GeV","mt2w_300_gev",100, 100.0, 490.0);
	TH1 *histMT2W_corte_500_GeV = new TH1F("MT2W_500_GeV","mt2w_500_gev",100, 100.0, 490.0);
	TH1 *histMT2BL = new TH1F("MT2BL","mt2bl",100, 100.0, 490.0);
	TH1 *histMT2BL_corte_100_GeV = new TH1F("MT2BL_100_GeV","mt2bl_100_gev",100, 100.0, 490.0);
	TH1 *histMT2BL_corte_300_GeV = new TH1F("MT2BL_300_GeV","mt2bl_300_gev",100, 100.0, 490.0);
	TH1 *histMT2BL_corte_500_GeV = new TH1F("MT2BL_500_GeV","mt2bl_500_gev",100, 100.0, 490.0);
	TH1 *h_MuonPT = new TH1F("MUON_PT","muon_pt", 20,0.0,200.0);

	//Correlations of ISR-JET with other observables

	
	TH1 *h_Delta_eta_bh_ISR = new TH1F("Delta_eta_bh_isr","Delta_Eta_bh_isr",500, -20.0, 20.0);
	TH1 *h_Delta_phi_bh_ISR = new TH1F("h_Delta_phi_bh_ISR","h_Delta_PHI_bh_ISR",500, -100.0, 20.0);
	TH1 *h_Delta_PT_bh_ISR = new TH1F("h_Delta_pt_bh_ISR","h_Delta_PT_bh_ISR",500, 0.0, 100.0);
	TH1 *h_Delta_eta_bl_ISR = new TH1F("Delta_eta_bl_isr","Delta_Eta_bl_isr",500, -20.0, 20.0);
	TH1 *h_Delta_phi_bl_ISR = new TH1F("h_Delta_phi_bl_ISR","h_Delta_PHI_bl_ISR",500, -5.0, 20.0);
	TH1 *h_Delta_PT_bl_ISR = new TH1F("h_Delta_pt_bl_ISR","h_Delta_PT_bl_ISR",500, 0.0, 100.0);
	TH1 *h_Cociente_PT_Muon_PT_ISR = new TH1F("Cociente_PT_Muon_PT_ISR","Cociente_pt_muon_pt_ISR",100, 0.0, 20.0);
	TH1 *h_Delta_Eta_Muon_ISR = new TH1F("Delta_Eta_Muon_ISR","Delta_eta_muon_ISR",100, -20.0, 20.0);
	TH1 *h_Delta_Phi_Muon_ISR = new TH1F("Delta_Phi_Muon_ISR","Delta_phi_muon_ISR",100, -20.0, 20.0);
	TH1 *h_Cociente_PT_Electron_PT_ISR = new TH1F("Cociente_PT_Electron_PT_ISR","Cociente_pt_electron_pt_ISR",100, 0.0, 20.0);
	TH1 *h_Delta_Eta_Electron_ISR = new TH1F("Delta_Eta_Electron_ISR","Delta_eta_electron_ISR",100, -12.0, 12.0);
	TH1 *h_Delta_Phi_Electron_ISR = new TH1F("Delta_Phi_Electron_ISR","Delta_phi_electron_ISR",100, -9.0, 9.0);
	TH1 *h_Delta_Eta_ISR_MissingET = new TH1F("Delta_Eta_ISR_MissingET","Delta_eta_ISR_missinget",100, -12.0, 12.0);
	TH1 *h_Delta_Phi_ISR_MissingET = new TH1F("Delta_Phi_ISR_MissingET","Delta_phi_ISR_missinget",100, -9.0, 9.0);
	TH1 *h_Cociente_MissingET_PT_ISR = new TH1F("Cociente_MissingET_PT_ISR","Cociente_missinget_pt_ISR",100, 0.0, 20.0);
        TH2 *h_PT_Muon_MissingET = new TH2F("PT_Muon_MissingEt","Missinget_muon_pt", 500, 0.0, 200, 500, 0.0, 100.0);
	TH2 *h_PT_Electron_MissingET = new TH2F("PT_Electron_MissingEt","Missinget_electron_pt", 500, 0.0, 200, 500, 0.0, 100.0);
        TH2 *h_PT_Muon_Delta_Eta_Muon_ISR = new TH2F("PT_Muon_Delta_Eta_Muon_ISR","pt_muon_delta_eta_muon_ISR", 500, 0.0, 200, 500, 0.0, 15.0);
	TH2 *h_PT_Electron_Delta_Eta_Electron_ISR = new TH2F("PT_Electron_Delta_Eta_Electron_ISR","pt_electron_delta_eta_electron_ISR", 500, 0.0, 200, 500, 0.0, 15.0);
	TH2 *h_Delta_Eta_Muon_ISR_Delta_Eta_ISR_MissingET = new TH2F("Delta_Eta_Muon_ISR_Delta_Eta_ISR_MissingET","delta_eta_muon_ISR_delta_eta_ISR_missinget", 500, 0.0, 200, 500, 0.0, 15.0);
	TH2 *h_Delta_Eta_Electron_ISR_Delta_Eta_ISR_MissingET = new TH2F("Delta_Eta_Electron_ISR_Delta_Eta_ISR_MissingET","delta_eta_electron_ISR_delta_eta_ISR_missinget", 500, 0.0, 200, 500, 0.0, 15.0);
	TH2 *h_PT_Muon_PT_ISR = new TH2F("PT_Muon_PT_ISR","pt_muon_pt_ISR", 500, 0.0, 200, 500, 0.0, 100.0);
	TH2 *h_PT_Electron_PT_ISR = new TH2F("PT_Electron_PT_ISR","pt_electron_pt_ISR", 500, 0.0, 200, 500, 0.0, 100.0);
	TH2 *h_RM_Delta_Phi_ISR_MissingET = new TH2F("RM_PHI_Delta_Phi_ISR_MissingET","rm_phi_delta_phi_isr_missinget", 100, 0.0, 5, 100, 0.0, 5.0);
	TH1 *h_RM = new TH1F("RM","rm",100, 0.0, 5.0);
	TH1 *h_RM_PHI = new TH1F("RM_PHI","rm_phi",100, 0.0, 5.0);

	//Correlations of ISR-JET with other observables with PT cuts in ISR
	//Corte de 100 GeV

	TH1 *h_Delta_eta_bh_ISR_corte_100_GeV = new TH1F("Delta_eta_bh_isr_corte_100_GeV","Delta_Eta_bh_isr_corte_100_GeV",100, -20.0, 20.0);
	TH1 *h_Delta_phi_bh_ISR_corte_100_GeV = new TH1F("h_Delta_phi_bh_ISR-corte_100_GeV","h_Delta_PHI_bh_ISR_corte_100_GeV",100, -5.0, 20.0);
	TH1 *h_Delta_PT_bh_ISR_corte_100_GeV = new TH1F("h_Delta_pt_bh_ISR_corte_100_GeV","h_Delta_PT_bh_ISR_corte_100_GeV",100, 0.0, 100.0);
	TH1 *h_Delta_eta_bl_ISR_corte_100_GeV = new TH1F("Delta_eta_bl_isr_corte_100_GeV","Delta_Eta_bl_isr_corte_100_GeV",100, -20.0, 20.0);
	TH1 *h_Delta_phi_bl_ISR_corte_100_GeV = new TH1F("h_Delta_phi_bl_ISR_corte_100_GeV","h_Delta_PHI_bl_ISR_corte_100_GeV",100, -5.0, 20.0);
	TH1 *h_Delta_PT_bl_ISR_corte_100_GeV = new TH1F("h_Delta_pt_bl_ISR_corte_100_GeV","h_Delta_PT_bl_ISR_corte_100_GeV",100, 0.0, 100.0);
	TH1 *h_Cociente_PT_Muon_PT_ISR_corte_100_GeV = new TH1F("Cociente_PT_Muon_PT_ISR_corte_100_GeV","Cociente_pt_muon_pt_ISR_corte_100_GeV",100, 0.0, 20.0);
	TH1 *h_Delta_Eta_Muon_ISR_corte_100_GeV = new TH1F("Delta_Eta_Muon_ISR_corte_100_GeV","Delta_eta_muon_ISR_corte_100_GeV",100, -20.0, 20.0);
	TH1 *h_Delta_Phi_Muon_ISR_corte_100_GeV = new TH1F("Delta_Phi_Muon_ISR_corte_100_GeV","Delta_phi_muon_ISR_corte_100_GeV",100, -20.0, 20.0);
	TH1 *h_Cociente_PT_Electron_PT_ISR_corte_100_GeV = new TH1F("Cociente_PT_Electron_PT_ISR_corte_100_GeV","Cociente_pt_electron_pt_ISR_corte_100_GeV",100, 0.0, 20.0);
	TH1 *h_Delta_Eta_Electron_ISR_corte_100_GeV = new TH1F("Delta_Eta_Electron_ISR_corte_100_GeV","Delta_eta_electron_ISR_corte_100_GeV",100, -20.0, 20.0);
	TH1 *h_Delta_Phi_Electron_ISR_corte_100_GeV = new TH1F("Delta_Phi_Electron_ISR_corte_100_GeV","Delta_phi_electron_ISR_corte_100_GeV",100, -20.0, 20.0);
	TH1 *h_Delta_Eta_ISR_MissingET_corte_100_GeV = new TH1F("Delta_Eta_ISR_MissingET_corte_100_GeV","Delta_eta_ISR_missinget_corte_100_GeV",100, -12.0, 12.0);
	TH1 *h_Delta_Phi_ISR_MissingET_corte_100_GeV = new TH1F("Delta_Phi_ISR_MissingET_corte_100_GeV","Delta_phi_ISR_missinget_corte_100_GeV",100, 0.0, 6.0);
	TH1 *h_Cociente_MissingET_PT_ISR_corte_100_GeV = new TH1F("Cociente_MissingET_PT_ISR_corte_100_GeV","Cociente_missinget_pt_ISR_corte_100_GeV",100, 0.0, 20.0);
        TH2 *h_PT_Muon_MissingET_corte_100_GeV = new TH2F("PT_Muon_MissingEt_corte_100_GeV","Missinget_muon_pt_corte_100_GeV", 500, 0.0, 200, 500, 0.0, 100.0);
	TH2 *h_PT_Electron_MissingET_corte_100_GeV = new TH2F("PT_Electron_MissingEt_corte_100_GeV","Missinget_electron_pt_corte_100_GeV", 500, 0.0, 200, 500, 0.0, 100.0);
        TH2 *h_PT_Muon_Delta_Eta_Muon_ISR_corte_100_GeV = new TH2F("PT_Muon_Delta_Eta_Muon_ISR_corte_100_GeV","pt_muon_delta_eta_muon_ISR_corte_100_GeV", 500, 0.0, 200, 500, 0.0, 15.0);
	TH2 *h_PT_Electron_Delta_Eta_Electron_ISR_100_GeV = new TH2F("PT_Electron_Delta_Eta_Electron_ISR_corte_100_GeV","pt_electron_delta_eta_electron_ISR_corte_100_GeV", 500, 0.0, 200, 500, 0.0, 15.0);
	TH2 *h_Delta_Eta_Muon_ISR_Delta_Eta_ISR_MissingET_corte_100_GeV = new TH2F("Delta_Eta_Muon_ISR_Delta_Eta_ISR_MissingET_corte_100_GeV","delta_eta_muon_ISR_delta_eta_ISR_missinget_corte_100_GeV", 500, 0.0, 200, 500, 0.0, 15.0);
	TH2 *h_Delta_Eta_Electron_ISR_Delta_Eta_ISR_MissingET_corte_100_GeV = new TH2F("Delta_Eta_Electron_ISR_Delta_Eta_ISR_MissingET_corte_100_GEV","delta_eta_electron_ISR_delta_eta_ISR_missinget_corte_100_GeV", 500, 0.0, 200, 500, 0.0, 15.0);
	TH2 *h_PT_Muon_PT_ISR_corte_100_GeV = new TH2F("PT_Muon_PT_ISR_corte_100_GeV","pt_muon_pt_ISR_corte_100_GeV", 500, 0.0, 200, 500, 0.0, 100.0);
	TH2 *h_PT_Electron_PT_ISR_corte_100_GeV = new TH2F("PT_Electron_PT_ISR_corte_100_GeV","pt_electron_pt_ISR_corte_100_GeV", 500, 0.0, 200, 500, 0.0, 100.0);
	TH2 *h_RM_Delta_Phi_ISR_MissingET_corte_100_GeV = new TH2F("RM_PHI_Delta_Phi_ISR_MissingET_corte_100_GeV","rm_phi_delta_phi_isr_missinget_corte_100_GeV", 100, 0.0, 5, 100, 0.0, 5.0);
	TH1 *h_RM_corte_100_GeV = new TH1F("RM_corte_100_GeV","rm_corte_100_GeV",100, 0.0, 5.0);
	TH1 *h_RM_PHI_corte_100_GeV = new TH1F("RM_PHI_corte_100_GeV","rm_phi_corte_100_GeV",100, 0.0, 5.0);


	//Corte de 300 GeV
	TH1 *h_Delta_eta_bh_ISR_corte_300_GeV = new TH1F("Delta_eta_bh_isr_corte_300_GeV","Delta_Eta_bh_isr_corte_300_GeV",100, -20.0, 20.0);
	TH1 *h_Delta_phi_bh_ISR_corte_300_GeV = new TH1F("h_Delta_phi_bh_ISR-corte_300_GeV","h_Delta_PHI_bh_ISR_corte_300_GeV",100, -5.0, 20.0);
	TH1 *h_Delta_PT_bh_ISR_corte_300_GeV = new TH1F("h_Delta_pt_bh_ISR_corte_300_GeV","h_Delta_PT_bh_ISR_corte_300_GeV",100, 0.0, 100.0);
	TH1 *h_Delta_eta_bl_ISR_corte_300_GeV = new TH1F("Delta_eta_bl_isr_corte_300_GeV","Delta_Eta_bl_isr_corte_300_GeV",100, -20.0, 20.0);
	TH1 *h_Delta_phi_bl_ISR_corte_300_GeV = new TH1F("h_Delta_phi_bl_ISR_corte_300_GeV","h_Delta_PHI_bl_ISR_corte_300_GeV",100, -5.0, 20.0);
	TH1 *h_Delta_PT_bl_ISR_corte_300_GeV = new TH1F("h_Delta_pt_bl_ISR_corte_300_GeV","h_Delta_PT_bl_ISR_corte_300_GeV",100, 0.0, 100.0);
	TH1 *h_Cociente_PT_Muon_PT_ISR_corte_300_GeV = new TH1F("Cociente_PT_Muon_PT_ISR_corte_300_GeV","Cociente_pt_muon_pt_ISR_corte_300_GeV",100, 0.0, 20.0);
	TH1 *h_Delta_Eta_Muon_ISR_corte_300_GeV = new TH1F("Delta_Eta_Muon_ISR_corte_300_GeV","Delta_eta_muon_ISR_corte_300_GeV",100, -20.0, 20.0);
	TH1 *h_Delta_Phi_Muon_ISR_corte_300_GeV = new TH1F("Delta_Phi_Muon_ISR_corte_300_GeV","Delta_phi_muon_ISR_corte_300_GeV",100, -20.0, 20.0);
	TH1 *h_Cociente_PT_Electron_PT_ISR_corte_300_GeV = new TH1F("Cociente_PT_Electron_PT_ISR_corte_300_GeV","Cociente_pt_electron_pt_ISR_corte_300_GeV",100, 0.0, 20.0);
	TH1 *h_Delta_Eta_Electron_ISR_corte_300_GeV = new TH1F("Delta_Eta_Electron_ISR_corte_300_GeV","Delta_eta_electron_ISR_corte_300_GeV",100, -20.0, 20.0);
	TH1 *h_Delta_Phi_Electron_ISR_corte_300_GeV = new TH1F("Delta_Phi_Electron_ISR_corte_300_GeV","Delta_phi_electron_ISR_corte_300_GeV",100, -20.0, 20.0);
	TH1 *h_Delta_Eta_ISR_MissingET_corte_300_GeV = new TH1F("Delta_Eta_ISR_MissingET_corte_300_GeV","Delta_eta_ISR_missinget_corte_300_GeV",100, -12.0, 12.0);
	TH1 *h_Delta_Phi_ISR_MissingET_corte_300_GeV = new TH1F("Delta_Phi_ISR_MissingET_corte_300_GeV","Delta_phi_ISR_missinget_corte_300_GeV",100, 0.0, 6.0);
	TH1 *h_Cociente_MissingET_PT_ISR_corte_300_GeV = new TH1F("Cociente_MissingET_PT_ISR_corte_300_GeV","Cociente_missinget_pt_ISR_corte_300_GeV",100, 0.0, 20.0);
        TH2 *h_PT_Muon_MissingET_corte_300_GeV = new TH2F("PT_Muon_MissingEt_corte_300_GeV","Missinget_muon_pt_corte_300_GeV", 500, 0.0, 200, 500, 0.0, 100.0);
	TH2 *h_PT_Electron_MissingET_corte_300_GeV = new TH2F("PT_Electron_MissingEt_corte_300_GeV","Missinget_electron_pt_corte_300_GeV", 500, 0.0, 200, 500, 0.0, 100.0);
        TH2 *h_PT_Muon_Delta_Eta_Muon_ISR_corte_300_GeV = new TH2F("PT_Muon_Delta_Eta_Muon_ISR_corte_300_GeV","pt_muon_delta_eta_muon_ISR_corte_300_GeV", 500, 0.0, 200, 500, 0.0, 15.0);
	TH2 *h_PT_Electron_Delta_Eta_Electron_ISR_300_GeV = new TH2F("PT_Electron_Delta_Eta_Electron_ISR_corte_300_GeV","pt_electron_delta_eta_electron_ISR_corte_300_GeV", 500, 0.0, 200, 500, 0.0, 15.0);
	TH2 *h_Delta_Eta_Muon_ISR_Delta_Eta_ISR_MissingET_corte_300_GeV = new TH2F("Delta_Eta_Muon_ISR_Delta_Eta_ISR_MissingET_corte_300_GeV","delta_eta_muon_ISR_delta_eta_ISR_missinget_corte_300_GeV", 500, 0.0, 200, 500, 0.0, 15.0);
	TH2 *h_Delta_Eta_Electron_ISR_Delta_Eta_ISR_MissingET_corte_300_GeV = new TH2F("Delta_Eta_Electron_ISR_Delta_Eta_ISR_MissingET_corte_300_GEV","delta_eta_electron_ISR_delta_eta_ISR_missinget_corte_300_GeV", 500, 0.0, 200, 500, 0.0, 15.0);
	TH2 *h_PT_Muon_PT_ISR_corte_300_GeV = new TH2F("PT_Muon_PT_ISR_corte_300_GeV","pt_muon_pt_ISR_corte_300_GeV", 500, 0.0, 200, 500, 0.0, 100.0);
	TH2 *h_PT_Electron_PT_ISR_corte_300_GeV = new TH2F("PT_Electron_PT_ISR_corte_300_GeV","pt_electron_pt_ISR_corte_300_GeV", 500, 0.0, 200, 500, 0.0, 100.0);
	TH2 *h_RM_Delta_Phi_ISR_MissingET_corte_300_GeV = new TH2F("RM_PHI_Delta_Phi_ISR_MissingET_corte_300_GeV","rm_phi_delta_phi_isr_missinget_corte_300_GeV", 100, 0.0, 5, 100, 0.0, 5.0);
	TH1 *h_RM_corte_300_GeV = new TH1F("RM_corte_300_GeV","rm_corte_300_GeV",100, 0.0, 5.0);
	TH1 *h_RM_PHI_corte_300_GeV = new TH1F("RM_PHI_corte_300_GeV","rm_phi_corte_300_GeV",100, 0.0, 5.0);


	//Corte de 500 GeV

	TH1 *h_Delta_eta_bh_ISR_corte_500_GeV = new TH1F("Delta_eta_bh_isr_corte_500_GeV","Delta_Eta_bh_isr_corte_500_GeV",100, -20.0, 20.0);
	TH1 *h_Delta_phi_bh_ISR_corte_500_GeV = new TH1F("h_Delta_phi_bh_ISR-corte_500_GeV","h_Delta_PHI_bh_ISR_corte_500_GeV",100, -5.0, 20.0);
	TH1 *h_Delta_PT_bh_ISR_corte_500_GeV = new TH1F("h_Delta_pt_bh_ISR_corte_500_GeV","h_Delta_PT_bh_ISR_corte_500_GeV",100, 0.0, 100.0);
	TH1 *h_Delta_eta_bl_ISR_corte_500_GeV = new TH1F("Delta_eta_bl_isr_corte_500_GeV","Delta_Eta_bl_isr_corte_500_GeV",100, -20.0, 20.0);
	TH1 *h_Delta_phi_bl_ISR_corte_500_GeV = new TH1F("h_Delta_phi_bl_ISR_corte_500_GeV","h_Delta_PHI_bl_ISR_corte_500_GeV",100, -5.0, 20.0);
	TH1 *h_Delta_PT_bl_ISR_corte_500_GeV = new TH1F("h_Delta_pt_bl_ISR_corte_500_GeV","h_Delta_PT_bl_ISR_corte_500_GeV",100, 0.0, 100.0);
	TH1 *h_Cociente_PT_Muon_PT_ISR_corte_500_GeV = new TH1F("Cociente_PT_Muon_PT_ISR_corte_500_GeV","Cociente_pt_muon_pt_ISR_corte_500_GeV",100, 0.0, 20.0);
	TH1 *h_Delta_Eta_Muon_ISR_corte_500_GeV = new TH1F("Delta_Eta_Muon_ISR_corte_500_GeV","Delta_eta_muon_ISR_corte_500_GeV",100, -20.0, 20.0);
	TH1 *h_Delta_Phi_Muon_ISR_corte_500_GeV = new TH1F("Delta_Phi_Muon_ISR_corte_500_GeV","Delta_phi_muon_ISR_corte_500_GeV",100, -20.0, 20.0);
	TH1 *h_Cociente_PT_Electron_PT_ISR_corte_500_GeV = new TH1F("Cociente_PT_Electron_PT_ISR_corte_500_GeV","Cociente_pt_electron_pt_ISR_corte_500_GeV",100, 0.0, 20.0);
	TH1 *h_Delta_Eta_Electron_ISR_corte_500_GeV = new TH1F("Delta_Eta_Electron_ISR_corte_500_GeV","Delta_eta_electron_ISR_corte_500_GeV",100, -20.0, 20.0);
	TH1 *h_Delta_Phi_Electron_ISR_corte_500_GeV = new TH1F("Delta_Phi_Electron_ISR_corte_500_GeV","Delta_phi_electron_ISR_corte_500_GeV",100, -20.0, 20.0);
	TH1 *h_Delta_Eta_ISR_MissingET_corte_500_GeV = new TH1F("Delta_Eta_ISR_MissingET_corte_500_GeV","Delta_eta_ISR_missinget_corte_500_GeV",100, -12.0, 12.0);
	TH1 *h_Delta_Phi_ISR_MissingET_corte_500_GeV = new TH1F("Delta_Phi_ISR_MissingET_corte_500_GeV","Delta_phi_ISR_missinget_corte_500_GeV",100, 0.0, 6.0);
	TH1 *h_Cociente_MissingET_PT_ISR_corte_500_GeV = new TH1F("Cociente_MissingET_PT_ISR_corte_500_GeV","Cociente_missinget_pt_ISR_corte_500_GeV",100, 0.0, 20.0);
        TH2 *h_PT_Muon_MissingET_corte_500_GeV = new TH2F("PT_Muon_MissingEt_corte_500_GeV","Missinget_muon_pt_corte_500_GeV", 500, 0.0, 200, 500, 0.0, 100.0);
	TH2 *h_PT_Electron_MissingET_corte_500_GeV = new TH2F("PT_Electron_MissingEt_corte_500_GeV","Missinget_electron_pt_corte_500_GeV", 500, 0.0, 200, 500, 0.0, 100.0);
        TH2 *h_PT_Muon_Delta_Eta_Muon_ISR_corte_500_GeV = new TH2F("PT_Muon_Delta_Eta_Muon_ISR_corte_500_GeV","pt_muon_delta_eta_muon_ISR_corte_500_GeV", 500, 0.0, 200, 500, 0.0, 15.0);
	TH2 *h_PT_Electron_Delta_Eta_Electron_ISR_500_GeV = new TH2F("PT_Electron_Delta_Eta_Electron_ISR_corte_500_GeV","pt_electron_delta_eta_electron_ISR_corte_500_GeV", 500, 0.0, 200, 500, 0.0, 15.0);
	TH2 *h_Delta_Eta_Muon_ISR_Delta_Eta_ISR_MissingET_corte_500_GeV = new TH2F("Delta_Eta_Muon_ISR_Delta_Eta_ISR_MissingET_corte_500_GeV","delta_eta_muon_ISR_delta_eta_ISR_missinget_corte_500_GeV", 500, 0.0, 200, 500, 0.0, 15.0);
	TH2 *h_Delta_Eta_Electron_ISR_Delta_Eta_ISR_MissingET_corte_500_GeV = new TH2F("Delta_Eta_Electron_ISR_Delta_Eta_ISR_MissingET_corte_500_GEV","delta_eta_electron_ISR_delta_eta_ISR_missinget_corte_500_GeV", 500, 0.0, 200, 500, 0.0, 15.0);
	TH2 *h_PT_Muon_PT_ISR_corte_500_GeV = new TH2F("PT_Muon_PT_ISR_corte_500_GeV","pt_muon_pt_ISR_corte_500_GeV", 500, 0.0, 200, 500, 0.0, 100.0);
	TH2 *h_PT_Electron_PT_ISR_corte_500_GeV = new TH2F("PT_Electron_PT_ISR_corte_500_GeV","pt_electron_pt_ISR_corte_500_GeV", 500, 0.0, 200, 500, 0.0, 100.0);

	TH2 *h_RM_Delta_Phi_ISR_MissingET_corte_500_GeV = new TH2F("RM_PHI_Delta_Phi_ISR_MissingET_corte_500_GeV","rm_phi_delta_phi_isr_missinget_corte_500_GeV", 100, 0.0, 5, 100, 0.0, 5.0);
	TH1 *h_RM_corte_500_GeV = new TH1F("RM_corte_500_GeV","rm_corte_500_GeV",100, 0.0, 5.0);
	TH1 *h_RM_PHI_corte_500_GeV = new TH1F("RM_PHI_corte_500_GeV","rm_phi_corte_500_GeV",100, 0.0, 5.0);

		
	//
	Int_t cut_cntr[13]={13*0}; // event counter after different cuts
	Int_t MET_cntr[40]={40*0}; // event counter after different cuts
        Int_t rmphi_cntr[30]={30*0};
        Int_t rm_cntr[30]={30*0};
        Int_t mt2_cntr[30]={30*0};

	//cut_cntr[0]= number of events read, cun_cntr[1]
	int n5jets=0;

	// Cycle over several runs . iRun corresponds to the seed of the current run

	//Stop 200 GeV y LSP 20 GeV
	for(int iRun = 1; iRun < 22; iRun ++){
	
	//Stop 200 GeV y LSP 35 GeV
	//for(int iRun = 1; iRun < 35; iRun ++){

	//Stop 300 GeV y LSP 135 GeV
	//for(int iRun = 1; iRun < 21; iRun ++){

	//Stop 300 GeV y LSP 120 GeV
	//for(int iRun = 1; iRun < 2; iRun ++){

	//Stop 400 GeV y LSP 220 GeV
	//for(int iRun = 1; iRun < 2; iRun ++){
	
	//Stop 400 GeV y LSP 235 GeV
	//for(int iRun = 1; iRun < 9; iRun ++){
	
	//Stop 500 GeV y LSP 320 GeV
	//for(int iRun = 1; iRun < 2; iRun ++){
	
	//Stop 400 GeV y LSP 100 GeV
	//for(int iRun = 1; iRun < 11; iRun ++){
	
	//Background semileptónico
	//for(int iRun = 1; iRun < 301; iRun ++){
	  
	//Background dileptónico
	//for(int iRun = 1; iRun < 97; iRun ++){
	  
	  
	// Create chains of root trees
	  TChain chain_Delphes("Delphes");
	  
	  // Loading simulations from Delphes
	  Char_t unidad = 0x30 + iRun%10;
	  Char_t decena = 0x30 + int(iRun/10)%10;
	  Char_t centena = 0x30 + int(iRun/100)%10;
	  
	  current_folder[current_folder.size()-4] = centena;
	  current_folder[current_folder.size()-3] = decena;
	  current_folder[current_folder.size()-2] = unidad;
	  matching_name[matching_name.size()-6] = centena;
	  matching_name[matching_name.size()-5] = decena;
	  matching_name[matching_name.size()-4] = unidad;
	  
	  string file_pythia_str = head_folder + current_folder + "Events/run_01/output_pythia8.root";
		
	  string file_delphes_str = head_folder + current_folder + "Events/run_01/output_delphes.root";
	  Char_t *file_delphes = (Char_t *) file_delphes_str.c_str();
	  
	  cout << "\nReading the file: \nDelphes: " << file_delphes << endl;
	  
	  chain_Delphes.Add(file_delphes);
	  // Objects of class ExRootTreeReader for reading the information
	  ExRootTreeReader *treeReader_Delphes = new ExRootTreeReader(&chain_Delphes);
	  
	  Long64_t numberOfEntries = treeReader_Delphes->GetEntries();
	  
		// Get pointers to branches used in this analysis
	  TClonesArray *branchJet = treeReader_Delphes->UseBranch("Jet");
	  TClonesArray *branchMissingET = treeReader_Delphes->UseBranch("MissingET");
	  TClonesArray *branchElectron = treeReader_Delphes->UseBranch("Electron");
	  TClonesArray *branchMuon = treeReader_Delphes->UseBranch("Muon");
	  
	  cout << endl;
	  cout << " Number of Entries Delphes = " << numberOfEntries << endl;
	  cout << endl;
	  
	  // particles, jets and vectors
	  MissingET *METpointer;
	  Muon *Muonpointer;
	  Electron *Electronpointer;
	  Jet *Jetpointer;
	  TLorentzVector *vect_currentJet = new TLorentzVector;
	  TLorentzVector *vect_auxJet = new TLorentzVector;
	  TLorentzVector *vect_leading = new TLorentzVector;
	  Jet *currentJet = new Jet;
	  Jet *auxJet = new Jet;
	  TRefArray array_temp;
	  
	  // Temporary variables
	  Double_t MET = 0.0; // Missing transverse energy
	  Double_t delta_phi = 0.0; // difference between the phi angle of MET and the jet
	  Double_t transverse_mass = 0.0; // Transverse mass
	  Double_t HT = 0.0; // Sum of jets' PT
	  Double_t HT_R1 = 0.0; // Sum of jets' PT which are in the same hemisphere of the ISR jet hemisphere
	  Double_t HT_R2 = 0.0; // Sum of jets' PT which are in the opposite hemisphere of the ISR jet hemisphere
	  Double_t ISR_Eta = 0.0; // Pseudorapidity of the ISR jet
	  Int_t number_Btags = 0; // Number of B jets per event
	  Int_t ISR_Btags = 0; // Number of BTags which are also ISR jets
	  Double_t delta_PT_jet = 0.0; // |PT-<PT>|
	  Double_t PT_sum = 0.0; // sum(PT)
	  Double_t PT_aver = 0.0; // <PT>
	  Double_t Delta_eta_aver = 0.0; // sum_i|eta-eta_i|/(Nj-1)
	  Double_t Delta_phi_sum = 0.0; // sum delta_phi
	  Double_t Delta_phi_other_jets = 0.0; // Average of delta phi of other jets
	  Double_t PT_ratio = 0.0; // PT/PT_others
	  Double_t Eta_ratio = 0.0; // Eta/Eta_others
	  Double_t Eta_sum = 0.0; // sum(Eta)
	  Double_t Delta_R = 0.0; // Size of the jet
	  Double_t Delta_phi_ratio = 0.0; // Delta_phi/Delta_phi_others
	  Double_t Delta_PT_leading = 0.0; // PT - PT_leading
	  Double_t Delta_Eta_leading = 0.0; // |Eta - Eta_leading|

	  /*
	   * Some variables used through the code
	   */
		Int_t ISR_jets[numberOfEntries];
		Int_t NumJets = 0;

		string fileName_str = head_folder_binary + matching_name;

		Char_t * fileName = (Char_t *) fileName_str.c_str();
		
		ifstream ifs(fileName,ios::in | ios::binary);

		for (Int_t j = 0; j<numberOfEntries; j++){
		  ifs.read((Char_t *) (ISR_jets+j),sizeof(Int_t));
		}
		ifs.close();
		
		// Jet with greatest PT
		Double_t PT_max = 0;
		Int_t posLeadingPT = -1;
		Int_t ISR_greatest_PT = 0;
		Double_t MT_leading_jet = 0.0; // Transverse mass
		bool ISR_found =false;
		int bjet_cntr=0;
		int iwjets=0;
		int jet_isr=999;
                //Correlacion de variables

                double Delta_eta_bl_ISR=0;
		double Delta_eta_bh_ISR=0;
                double Delta_phi_bl_ISR=0;
		double Delta_phi_bh_ISR=0;
		double Delta_PT_bl_ISR=0;
                double Delta_PT_bh_ISR=0;      
                double rm=0;
                double rmphi=0;


		//cout <<"paso los tclones1"<< endl; 
		/*
		 * Main cycle of the program
		 */
		numberOfEntries = 100000;
		//numberOfEntries = 8400;
		for (Int_t entry = 0; entry < numberOfEntries; ++entry){
		  // Progress
		  if(numberOfEntries>10 && (entry%((int)numberOfEntries/10))==0.0){
		    cout<<"progress = "<<(entry*100/numberOfEntries)<<"%\t";
		    cout<< "Time :"<< (clock()-initialTime)/double_t(CLOCKS_PER_SEC)<<"s"<<endl;
		  }
		  //cout <<"paso los tclones2"<< endl;
		  // Load selected branches with data from specified event
		  treeReader_Delphes->ReadEntry(entry);
		  
		  ISR_found=false;
		  bjet_cntr=0;
		  iwjets=0;
		  jet_isr=0;
		  cut_cntr[0]++;//Numero de eventos procesados
		  // MET
		  METpointer = (MissingET*) branchMissingET->At(0);
		  MET = METpointer->MET;
		  h_MET->Fill(MET);
			
		  Electronpointer = (Electron*)branchElectron->At(0);
		  Muonpointer = (Muon*)branchMuon->At(0);
		  Jetpointer = (Jet*)branchJet->At(0);
		  
		  if(branchMuon->GetEntries() >0 )
		    {
		      
		      h_MuonPT->Fill(Muonpointer->PT);
		    }
		  NumJets=branchJet->GetEntries();
		  h_numberJet->Fill(NumJets);
		  
		  if (NumJets==5) n5jets++;
		  
		  // checking the ISR
		  if (ISR_jets[entry] == -1 || NumJets < 3)
		    continue;
		  
		  PT_max = 0;
		  posLeadingPT = -1;
		  HT = 0;
		  HT_R1 = 0;
		  HT_R2 = 0;
		  number_Btags = 0;
		  
		  delta_PT_jet = 0.0;
		  PT_aver = 0.0;
		  PT_sum = 0.0;
		  Delta_eta_aver = 0.0;
		  Delta_phi_sum = 0.0;
		  Delta_phi_other_jets = 0.0;
		  Delta_phi_ratio = 0.0;
		  Delta_PT_leading = 0.0;
		  Delta_Eta_leading = 0.0;

		  PT_ratio = 0.0;
		  Eta_ratio = 0.0;
		  Eta_sum = 0.0;
		  
		  Delta_R = 0.0;

		  if (ISR_jets[entry] >= NumJets){
		    cout << "Error en el matching" << endl;
		    return 1;
		  }
		  
		  // Preliminary for. It is used to calculate PT_aver and Delta_phi_sum
		  for (int iJet = 0; iJet<NumJets; iJet++){
		    currentJet = (Jet*) branchJet->At(iJet);
		    vect_currentJet->SetPtEtaPhiM(currentJet->PT,currentJet->Eta,currentJet->Phi,currentJet->Mass);
		    delta_phi = deltaAng(vect_currentJet->Phi(), METpointer->Phi);
		    PT_sum += vect_currentJet->Pt();
		    Eta_sum += vect_currentJet->Eta();
		    Delta_phi_sum += delta_phi;
		    // HT
		    HT += vect_currentJet->Pt();
		    // HT ratios
		    if((vect_currentJet->Eta()*ISR_Eta) > 0)
		      HT_R1 += vect_currentJet->Pt();
		    else
		      HT_R2 += vect_currentJet->Pt();
		    // PT Leading jet
		    if(PT_max < vect_currentJet->Pt()){
		      PT_max = vect_currentJet->Pt();
		      posLeadingPT = iJet;
		    }
		  }
		  
		  //PT_aver
		  PT_aver = PT_sum/NumJets;
		  
		  // Leading PT
		  currentJet = (Jet*) branchJet->At(posLeadingPT);
		  vect_leading->SetPtEtaPhiM(currentJet->PT,currentJet->Eta,currentJet->Phi,currentJet->Mass);
		  
		  // ISR jet
		  currentJet = (Jet*) branchJet->At(ISR_jets[entry]);
		  vect_currentJet->SetPtEtaPhiM(currentJet->PT,currentJet->Eta,currentJet->Phi,currentJet->Mass);
		  ISR_Eta = vect_currentJet->Eta();
      
		  int bjet[5]={5*0};
		  int wjets[5]={5*0}; 
		  
		  for (int iJet = 0; iJet<NumJets; iJet++){
		    currentJet = (Jet*) branchJet->At(iJet);
		    vect_currentJet->SetPtEtaPhiM(currentJet->PT,currentJet->Eta,currentJet->Phi,currentJet->Mass);
		    delta_phi = deltaAng(vect_currentJet->Phi(), METpointer->Phi);
		    //transverse_mass = sqrt(2*vect_currentJet->Pt()*MET*(1-cos(delta_phi)));
		    
		    // Correlated variables
		    delta_PT_jet = TMath::Abs(vect_currentJet->Pt()-PT_aver);
		    Delta_phi_other_jets = (Delta_phi_sum-delta_phi)/(NumJets-1);
		    PT_ratio = vect_currentJet->Pt()*(NumJets-1)/(PT_sum-vect_currentJet->Pt());
		    Eta_ratio = vect_currentJet->Eta()*(NumJets-1)/(Eta_sum-vect_currentJet->Eta());
		    Delta_phi_ratio = delta_phi*(NumJets-1)/(Delta_phi_sum-delta_phi);
		    
		    Delta_Eta_leading = TMath::Abs(vect_currentJet->Eta()-vect_leading->Eta());
		    Delta_PT_leading = vect_leading->Pt()-vect_currentJet->Pt();
		    
		    Delta_eta_aver = 0.0;
		    // For cycle used to calculate Delta_eta_aver
		    for(Int_t iJet2 = 0; iJet2<NumJets; iJet2++){
		      auxJet = (Jet*) branchJet->At(iJet2);
		      vect_auxJet->SetPtEtaPhiM(auxJet->PT,auxJet->Eta,auxJet->Phi,auxJet->Mass);
		      if (iJet2 != iJet) Delta_eta_aver += TMath::Abs(vect_auxJet->Eta()-vect_currentJet->Eta());
		    }
		    Delta_eta_aver = Delta_eta_aver/(NumJets-1);
		    Delta_R = sqrt(pow(currentJet->DeltaEta,2)+pow(currentJet->DeltaPhi,2));
		    
		    // Multiplicity
		    array_temp = (TRefArray) currentJet->Constituents;
		    
		    if (iJet != ISR_jets[entry]){ // Non ISR
		      h_jet_PT->Fill(vect_currentJet->Pt());
		      h_jet_Eta->Fill(vect_currentJet->Eta());
		      h_jet_Phi->Fill(vect_currentJet->Phi());
		      h_jet_DPhi_MET->Fill(delta_phi);
		      h_jet_MT->Fill(transverse_mass);
		      h_jet_Delta_PT->Fill(delta_PT_jet);
		      h_jet_Delta_Eta->Fill(Delta_eta_aver);
		      h_jet_DPhi_MET_other->Fill(Delta_phi_other_jets);
		      h_jet_PT_HT->Fill(vect_currentJet->Pt()/HT);
		      h_jet_multiplicity->Fill(array_temp.GetEntries());
		      h_jet_PT_over_PT_others->Fill(PT_ratio);
		      h_jet_Eta_over_Eta_others->Fill(Eta_ratio);
		      h_jet_DeltaR->Fill(Delta_R);
		      h_jet_DPhi_over_Phi_others->Fill(Delta_phi_ratio);
		      h_jet_Delta_PT_leading->Fill(Delta_PT_leading);
		      h_jet_Delta_Eta_leading->Fill(Delta_Eta_leading);
		      if (vect_currentJet->Pt()>240)
			h_jet_DPhi_MET_hpt->Fill(delta_phi);
		      h2_jet_PTEta->Fill(vect_currentJet->Pt(),vect_currentJet->Eta());
		      
		      // For testing creating histo
		      hist_jet_PT->Fill(vect_currentJet->Pt());
		      hist_jet_Abs_Eta->Fill(TMath::Abs(vect_currentJet->Eta()));
		      hist_jet_DPhi_MET->Fill(delta_phi);
		      hist_jet_PT_ratio->Fill(PT_ratio);
		      hist_jet_Delta_Eta->Fill(Delta_eta_aver);
		      hist_jet_DPhi_MET_other->Fill(Delta_phi_other_jets);
		      hist_jet_Delta_PT_leading->Fill(Delta_PT_leading);
		      hist_jet_Delta_Eta_leading->Fill(Delta_Eta_leading);
		      
		    }
		    
		    else{ //ISR
		      //cout <<"paso los tclones1"<< endl; 
		      
		      
		      //Masa invariante de los jets
		      
		      ISR_found=true;
		      jet_isr=iJet;
		      
		      h_ISR_PT->Fill(vect_currentJet->Pt());
		      h_ISR_Eta->Fill(vect_currentJet->Eta());
		      h_ISR_Phi->Fill(vect_currentJet->Phi());
		      h_ISR_DPhi_MET->Fill(delta_phi);
		      h_ISR_Eta_comp->Fill(vect_currentJet->Eta());
		      h_ISR_PT_comp->Fill(vect_currentJet->Pt());
		      h_ISR_DPhi_MET_comp->Fill(delta_phi);
		      h_ISR_Delta_PT->Fill(delta_PT_jet);
		      h_ISR_Delta_Eta->Fill(Delta_eta_aver);
		      h_ISR_DPhi_MET_other->Fill(Delta_phi_other_jets);
		      h_ISR_PT_HT->Fill(vect_currentJet->Pt()/HT);
		      h_ISR_multiplicity->Fill(array_temp.GetEntries());
		      h_ISR_PT_over_PT_others->Fill(PT_ratio);
		      h_ISR_Eta_over_Eta_others->Fill(Eta_ratio);
		      h_ISR_DeltaR->Fill(Delta_R);
		      h_ISR_DPhi_over_Phi_others->Fill(Delta_phi_ratio);
		      h_ISR_Delta_PT_leading->Fill(Delta_PT_leading);
		      h_ISR_Delta_Eta_leading->Fill(Delta_Eta_leading);
		      
		      double MET_sqrt_HT=transverse_mass/sqrt(HT);
		      h_ISR_MET_HT->Fill(MET_sqrt_HT);
		      if (vect_currentJet->Pt()>100)
			h_ISR_MET_HT1->Fill(MET_sqrt_HT);
		      if (vect_currentJet->Pt()>200)
			h_ISR_MET_HT2->Fill(MET_sqrt_HT);
		      if (vect_currentJet->Pt()>300){
			h_ISR_MET_HT3->Fill(MET_sqrt_HT);
			
		      }
			    
		      // For testing creating histo
		      hist_ISR_PT->Fill(vect_currentJet->Pt());
		      hist_ISR_Abs_Eta->Fill(TMath::Abs(vect_currentJet->Eta()));
		      hist_ISR_DPhi_MET->Fill(delta_phi);
		      hist_ISR_PT_ratio->Fill(PT_ratio);
		      hist_ISR_Delta_Eta->Fill(Delta_eta_aver);
		      hist_ISR_DPhi_MET_other->Fill(Delta_phi_other_jets);
		      hist_ISR_Delta_PT_leading->Fill(Delta_PT_leading);
		      hist_ISR_Delta_Eta_leading->Fill(Delta_Eta_leading);
		    } //else ISR
		    
			  // BTag
		    h_BTag->Fill(currentJet->BTag);
		    if (currentJet->BTag == 1){ // The current jet is B Tagged
		      h_BTag_PT->Fill(vect_currentJet->Pt());
		      h_BTag_Eta->Fill(vect_currentJet->Eta());
		      h_BTag_DPhi_MET->Fill(delta_phi);
		      number_Btags++;
		      
		      if (iJet == ISR_jets[entry]){ // If the ISR jet is also a B jet
			ISR_Btags++;
		      }
		      
		    }
		    
		    // Checking type of jets
		    if (currentJet->BTag>0){
		      bjet[bjet_cntr]= iJet;
		      bjet_cntr++;
		    }
		    else{ if (iJet!=jet_isr){
			wjets[iwjets]=iJet;
			iwjets++;
		      }
		    }
		    
			} // end loop jets
		  
		  
		  
		  bool  esolo=false;
		  bool musolo=false;
		  double PT_Muon;
		  double PT_Electron;
		  double PT_ISR;
		  double MissingET;
		  double transverse_mass_Muon=0;
		  double transverse_mass_Electron=0;
		  double transverse_mass_square_Muon;
		  double transverse_mass_square_Electron;
		  double Cociente_PT_Muon_PT_ISR;
		  double Delta_Eta_Muon_ISR;
		  double Delta_Phi_Muon_ISR;
		  double Cociente_PT_Electron_PT_ISR;
		  double Delta_Eta_Electron_ISR;
		  double Delta_Phi_Electron_ISR;
		  double Cociente_MissingET_PT_ISR;
		  double Delta_Eta_ISR_MissingET;
		  double Delta_Phi_ISR_MissingET;
		  Jet *isrjet = new Jet;
		  Jet *wjet1 = new Jet;
		  Jet *wjet2 = new Jet;
		  Jet *bjet1 = new Jet;
		  Jet *bjet2 = new Jet;
		  
		  

		  if ((branchElectron->GetEntries() == 1) && (branchMuon->GetEntries() == 0)) esolo=true;
		  if ((branchElectron->GetEntries() == 0) && (branchMuon->GetEntries() == 1)) musolo=true;

		  if (esolo||musolo) cut_cntr[1]++; //Numero de eventos con muon o electron
		  //if (ISR_found) cut_cntr[2]++;}
		  //if ((esolo||musolo)&&(bjet_cntr==1)) cut_cntr[2]++;
		  //if ((esolo||musolo)&&(bjet_cntr==1)&&(NumJets>3)) cut_cntr[3]++;
		  
		  if (ISR_found&&(esolo||musolo)){
		    
		    
		    double pl[4]; // lepton cuadrimomentum
		    double pb1[4]; // bjet1 cuadrimomentum, bjet1 belong to the leptonic branch
		    double pb2[4]; // bjet2 cuadrimomentum
		    isrjet = (Jet*) branchJet->At(jet_isr);
		    
		    cut_cntr[2]++;//Numero de eventos con ISR + lepton solo 
		    if (bjet_cntr==1){ cut_cntr[3]++;} //Numero de eventos con 1 bjet
			
		    if (bjet_cntr==2){	
		      cut_cntr[4]++;//Numero de eventos con dos bjets
		      bjet1 = (Jet*) branchJet->At(bjet[0]);
		      bjet2 = (Jet*) branchJet->At(bjet[1]);

		      if(iwjets==2){
		      cut_cntr[5]++;//Numero de eventos con dos wjets
		      wjet1 = (Jet*) branchJet->At(wjets[0]);
		      wjet2 = (Jet*) branchJet->At(wjets[1]);
	
		      if ((bjet1->PT>30)&&(bjet2->PT>30)){
		      cut_cntr[6]++;//Numero de eventos con 2 bjet y con corte en el PT > 30 GeV

		      		      
		      if(musolo)
			{
			  //Cociente de PT del Muón / PT del ISR 
			  Cociente_PT_Muon_PT_ISR = Muonpointer->PT/isrjet->PT;
			  h_Cociente_PT_Muon_PT_ISR->Fill(Cociente_PT_Muon_PT_ISR);
			  
			  
			  //Diferencia en ETA del Muón y del ISR
			  Delta_Eta_Muon_ISR = Muonpointer->Eta-isrjet->Eta;
			  h_Delta_Eta_Muon_ISR->Fill(Delta_Eta_Muon_ISR);
			  
			  
			  //Diferencia en PHI del Muón y del ISR
			  Delta_Phi_Muon_ISR = Muonpointer->Phi-isrjet->Phi;
			  h_Delta_Phi_Muon_ISR->Fill(Delta_Phi_Muon_ISR);
			  
			  //Histogramas en dos dimensiones del PT del muón y el MissingET
			  PT_Muon = Muonpointer->PT;
			  MissingET = METpointer->MET; 
			  h_PT_Muon_MissingET->Fill(PT_Muon,MissingET);

			  //Histogramas en dos dimensiones del PT del muón y la diferencia en ETA del muón y del ISR 
			  h_PT_Muon_Delta_Eta_Muon_ISR->Fill(PT_Muon,Delta_Eta_Muon_ISR);

			  //Histograma en dos dimensiones del PT del muón y el PT del ISR
			  PT_ISR = vect_currentJet->Pt();
			  h_PT_Muon_PT_ISR->Fill(PT_Muon,PT_ISR);
			  
			  //Masa Transversa
			  transverse_mass_Muon = TMath::  Sqrt ( Muonpointer->PT * METpointer->MET * (1 - cos(Muonpointer->Phi - METpointer->Phi)));
			  h_MT_Muon->Fill(transverse_mass_Muon);
			  //Masa transversa al cuadrado
			  transverse_mass_square_Muon = Muonpointer->PT * METpointer->MET * (1 - cos(Muonpointer->Phi - METpointer->Phi));
			  h_MT_square_Muon->Fill(transverse_mass_square_Muon);
			}
			    
		      if(esolo)
			{
			  //Cociente de PT del electrón / PT del ISR 
			  Cociente_PT_Electron_PT_ISR = Electronpointer->PT/isrjet->PT;
			  h_Cociente_PT_Electron_PT_ISR->Fill(Cociente_PT_Electron_PT_ISR);
				

			  //Diferencia en ETA del electrón y del ISR
			  Delta_Eta_Electron_ISR = Electronpointer->Eta-isrjet->Eta;
			  h_Delta_Eta_Electron_ISR->Fill(Delta_Eta_Electron_ISR);

			  //Diferencia en PHI del electrón y del ISR
			  Delta_Phi_Electron_ISR = Electronpointer->Phi-isrjet->Phi;
			  h_Delta_Phi_Electron_ISR->Fill(Delta_Phi_Electron_ISR);
			  
			  //Histogramas en dos dimensiones del PT del electrón y el MissingET
			  PT_Electron = Electronpointer->PT;
			  MissingET = METpointer->MET; 
			  h_PT_Electron_MissingET->Fill(PT_Electron,MissingET);
			  
			  //Histogramas en dos dimensiones del PT del electrón y la diferencia en ETA del electrón y del ISR 
			  h_PT_Electron_Delta_Eta_Electron_ISR->Fill(PT_Electron,Delta_Eta_Electron_ISR);

			  //Histograma en dos dimensiones del PT del electrón y el PT del ISR
			  PT_ISR = isrjet->PT;
			  h_PT_Electron_PT_ISR->Fill(PT_Electron,PT_ISR);
			  
			  //Masa transversa
			  transverse_mass_Electron = TMath::  Sqrt ( Electronpointer->PT * METpointer->MET * (1 - cos(Electronpointer->Phi - METpointer->Phi)));
			  h_MT_Electron->Fill(transverse_mass_Electron);
			  //Masa transversa al cuadrado
			  transverse_mass_square_Electron = Electronpointer->PT * METpointer->MET * (1 - cos(Electronpointer->Phi - METpointer->Phi));
			  h_MT_square_Electron->Fill(transverse_mass_square_Electron);	
			}
		    
		      //Cociente de la energia transversa faltante / PT del ISR
		      Cociente_MissingET_PT_ISR = METpointer->MET/isrjet->PT;
		      h_Cociente_MissingET_PT_ISR->Fill(Cociente_MissingET_PT_ISR);
		      
		      //Diferencia en ETA del ISR y del MissingET
		      Delta_Eta_ISR_MissingET = isrjet->Eta-METpointer->Eta;
		      h_Delta_Eta_ISR_MissingET->Fill(Delta_Eta_ISR_MissingET);
		      
		      //Diferencia en PHI del ISR y del MissingET
		      Delta_Phi_ISR_MissingET = TMath::Sqrt((isrjet->Phi-METpointer->Phi)*(isrjet->Phi-METpointer->Phi));
		      h_Delta_Phi_ISR_MissingET->Fill(Delta_Phi_ISR_MissingET);
		      
		      // Histograma en dos dimensiones del Delta Eta del leptón y del ISR en función del Delta Eta del ISR y del MissingET
		      h_Delta_Eta_Muon_ISR_Delta_Eta_ISR_MissingET->Fill(Delta_Eta_Muon_ISR,Delta_Eta_ISR_MissingET);
		      
		      TLorentzVector jetvec1, jetvec2, jetvec3;
		      TLorentzVector jetvecb1, jetvecb2, jetvecb3;
		      	     
		      
		      //cout<<"wjets[0]="<<wjets[0]<<", wjets[1]="<<wjets[1]<<endl;
		      //cout<<"bjet[0]="<<bjet[0]<<", bjet[1]="<<bjet[1]<<"; isr_jet="<<jet_isr<<endl;
	      
		      // Calculate the invariant mass of w jets
		      
		      jetvec1.SetPtEtaPhiM(wjet1->PT, wjet1->Eta, wjet1->Phi, 0.0);
		      jetvec2.SetPtEtaPhiM(wjet2->PT, wjet2->Eta, wjet2->Phi, 0.0);
		      
		      // Calculate the invariant mass of t jets
		      jetvecb1.SetPtEtaPhiM(bjet1->PT, bjet1->Eta, bjet1->Phi, 0.0);
		      jetvecb2.SetPtEtaPhiM(bjet2->PT, bjet2->Eta, bjet2->Phi, 0.0);
		      
		      double theta_bjet1 = 2*TMath::ATan(TMath::Exp(-bjet1->Eta));
		      double theta_bjet2 = 2*TMath::ATan(TMath::Exp(-bjet2->Eta));
		      double theta_wjet1 = 2*TMath::ATan(TMath::Exp(-wjet1->Eta));
		      double theta_wjet2 = 2*TMath::ATan(TMath::Exp(-wjet2->Eta));
		      
		      double cos_theta_bjet1 = TMath::Cos(theta_bjet1);
		      double cos_theta_bjet2 = TMath::Cos(theta_bjet2);
		      double cos_theta_wjet1 = TMath::Cos(theta_wjet1);
		      double cos_theta_wjet2 = TMath::Cos(theta_wjet2);

		      double sin_theta_bjet1 = TMath::Sin(theta_bjet1);
		      double sin_theta_bjet2 = TMath::Sin(theta_bjet2);
		      double sin_theta_wjet1 = TMath::Sin(theta_wjet1);
		      double sin_theta_wjet2 = TMath::Sin(theta_wjet2);
			    
		      double px_top1, px_top2;
		      double py_top1, py_top2;
		      double px_lepton, py_lepton, pz_lepton, E_lepton={4*0.0};
		      double theta_lepton;
		      double m_electron=0.511e-3; 
		      double m_muon=105e-3;
		      
		      double  pb1x,pb1y,pb1z,pb2x,pb2y,pb2z;
		      double  pwj1x,pwj1y,pwj1z,pwj2x,pwj2y,pwj2z;
		      double eb1, eb2,ewj1,ewj2;
		      double Dchi1, Dchi2;
		      double ejet1, pjet1,ejet2,pjet2;
		      
		      double invm1;
		      double invm2;
		      
		      // Dchi1=(jetvec1 + jetvec2+ jetvecb1).M();
		      //Dchi2=(jetvec1 + jetvec2+ jetvecb2).M();
		      
		      pb1x=bjet1->PT*TMath::Cos(bjet1->Phi);
		      pb1y=bjet1->PT*TMath::Sin(bjet1->Phi);	
		      pb1z=bjet1->PT/TMath::Sin(theta_bjet1);
		      eb1=TMath::Sqrt(TMath::Power(bjet1->Mass,2)+TMath::Power(pb1x,2)+TMath::Power(pb1y,2)+TMath::Power(pb1z,2));
		      
		      pb2x=bjet2->PT*TMath::Cos(bjet2->Phi);
		      pb2y=bjet2->PT*TMath::Sin(bjet2->Phi);	
		      pb2z=bjet2->PT/TMath::Sin(theta_bjet2);
		      eb2=TMath::Sqrt(TMath::Power(bjet2->Mass,2)+TMath::Power(pb2x,2)+TMath::Power(pb2y,2)+TMath::Power(pb2z,2));
		      
		      pwj1x=wjet1->PT*TMath::Cos(wjet1->Phi);
		      pwj1y=wjet1->PT*TMath::Sin(wjet1->Phi);
		      pwj1z=wjet1->PT/TMath::Sin(theta_wjet1);
		      ewj1=TMath::Sqrt(TMath::Power(wjet1->Mass,2)+TMath::Power(pwj1x,2)+TMath::Power(pwj1y,2)+TMath::Power(pwj1z,2));
		      
		      pwj2x=wjet2->PT*TMath::Cos(wjet2->Phi);
		      pwj2y=wjet2->PT*TMath::Sin(wjet2->Phi);	   
		      pwj2z=wjet2->PT/TMath::Sin(theta_wjet2);
		      ewj2=TMath::Sqrt(TMath::Power(wjet2->Mass,2)+TMath::Power(pwj2x,2)+TMath::Power(pwj2y,2)+TMath::Power(pwj2z,2));
		      
		      // Comprobación alternativa de la masa invariante usando Tlorentz Vector

		      invm1=(jetvec1 + jetvec2+ jetvecb1).M();
		      invm2=(jetvec1 + jetvec2+ jetvecb2).M();
			    
		      ejet1=(eb1+ewj1+ewj2)*(eb1+ewj1+ewj2);
		      pjet1=(pb1x+pwj1x+pwj2x)*(pb1x+pwj1x+pwj2x)+(pb1y+pwj1y+pwj2y)*(pb1y+pwj1y+pwj2y)+(pb1z+pwj1z+pwj2z)*(pb1z+pwj1z+pwj2z);
		      ejet2=(eb2+ewj1+ewj2)*(eb2+ewj1+ewj2);
		      pjet2=(pb2x+pwj1x+pwj2x)*(pb2x+pwj1x+pwj2x)+(pb2y+pwj1y+pwj2y)*(pb2y+pwj1y+pwj2y)+(pb2z+pwj1z+pwj2z)*(pb2z+pwj1z+pwj2z);      
		      Dchi1=sqrt(ejet1-pjet1);
		      Dchi2=sqrt(ejet2-pjet2);
		      /*
			cout<<" PT_bjet1="<<bjet1->PT<<" , eta_bjet1="<<bjet1->Eta<<" , phi_bjet1 ="<<bjet1->Phi<<" , mass_bjet1 ="<<bjet1->Mass<<endl;
			cout<<" Px_bjet1="<<pb1x<<" , Py_bjet1="<<pb1y<<" , Pz_bjet1 ="<<pb1z<<" , Theta_bjet1="<<theta_bjet1<<endl;
			cout<<" PT_bjet2="<<bjet2->PT<<" , eta_bjet2="<<bjet2->Eta<<" , phi_bjet2 ="<<bjet2->Phi<<" , mass_bjet2 ="<<bjet2->Mass<<endl;
			cout<<" Px_bjet2="<<pb2x<<" , Py_bjet2="<<pb2y<<" , Pz_bjet2 ="<<pb2z<<" , Theta_bjet2="<<theta_bjet2<<endl;
			cout<<" PT_wjet1="<<wjet1->PT<<" , eta_wjet1="<<wjet1->Eta<<" , phi_wjet1 ="<<wjet1->Phi<<" , mass_wjet1 ="<<wjet1->Mass<<endl;
			cout<<" Px_wjet1="<<pwj1x<<" , Py_wjet1="<<pwj1y<<" , Pz_wjet1 ="<<pwj1z<<" , Theta_wjet1="<<theta_wjet1<<endl; 
			cout<<" PT_wjet2="<<wjet2->PT<<" , eta_wjet2="<<wjet2->Eta<<" , phi_wjet2 ="<<wjet2->Phi<<" , mass_wjet2 ="<<wjet2->Mass<<endl;
			cout<<" Px_wjet2="<<pwj2x<<" , Py_wjet2="<<pwj2y<<" , Pz_wjet2 ="<<pwj2z<<" , Theta_wjet2="<<theta_wjet2<<endl; 
			cout<<" invm1="<<invm1<<" , invm2="<<invm2<<endl;
			cout<<"Dchi1="<<Dchi1<<" , Dchi2="<<Dchi2<<endl;
                        
		      */

		      //cout<<"invm1-dchi1="<<invm1-Dchi1<<"; invm2-dchi2="<<invm2-Dchi2<<endl;
		      
		      if (esolo) {
			theta_lepton=2*TMath::ATan(TMath::Exp(-Electronpointer->Eta));
			px_lepton=Electronpointer->PT*cos(Electronpointer->Phi);
			py_lepton=Electronpointer->PT*sin(Electronpointer->Phi);
			pz_lepton=Electronpointer->PT/sin(theta_lepton);
			E_lepton=sqrt(Electronpointer->PT*Electronpointer->PT + pz_lepton*pz_lepton + m_electron*m_electron);
			
		      }
		      
			    
		      if (musolo) {

			theta_lepton=2*TMath::ATan(TMath::Exp(-Muonpointer->Eta));
			px_lepton=Muonpointer->PT*cos(Muonpointer->Phi);
			py_lepton=Muonpointer->PT*sin(Muonpointer->Phi);
			pz_lepton=Muonpointer->PT/sin(theta_lepton);
			E_lepton=sqrt(Muonpointer->PT*Muonpointer->PT + pz_lepton*pz_lepton + m_muon*m_muon);
			      
		      }
		      
		      pl[0]=E_lepton;
		      pl[1]=px_lepton;
		      pl[2]=py_lepton;
		      pl[3]=pz_lepton;
		      if ((Dchi2-mtop)*(Dchi2-mtop)<(Dchi1-mtop)*(Dchi1-mtop)){
			
			//bjet2 is associated to hadronic w
			px_top1=bjet2->PT*cos(bjet2->Phi)+(wjet2->PT*cos(wjet2->Phi))+(wjet1->PT*cos(wjet1->Phi));
			py_top1=bjet2->PT*sin(bjet2->Phi)+(wjet2->PT*sin(wjet2->Phi))+(wjet1->PT*sin(wjet1->Phi));
			
			px_top2=bjet1->PT*cos(bjet1->Phi)+px_lepton;
			py_top2 =bjet1->PT*sin(bjet1->Phi)+py_lepton;
			
			pb1[1]=bjet1->PT*cos(bjet1->Phi);			      
			pb1[2]=bjet1->PT*sin(bjet1->Phi);				
			pb1[3]=bjet1->PT/sin(theta_bjet1);
			pb1[0]=sqrt(bjet1->PT*bjet1->PT+pb1[3]*pb1[3]+bjet1->Mass);
			
			pb2[1]=bjet2->PT*cos(bjet2->Phi);			      
			pb2[2]=bjet2->PT*sin(bjet2->Phi);				
			pb2[3]=bjet2->PT/sin(theta_bjet2);
			pb2[0]=sqrt(bjet2->PT*bjet2->PT+pb2[3]*pb2[3]+bjet2->Mass);
			//Rama hadronica
			Delta_eta_bh_ISR = bjet2->Eta-isrjet->Eta;
			Delta_phi_bh_ISR = bjet2->Phi-isrjet->Phi;
			Delta_PT_bh_ISR = bjet2->PT-isrjet->PT;
			//Rama leptonica
			Delta_eta_bl_ISR = bjet1->Eta-isrjet->Eta;
			Delta_phi_bl_ISR = bjet1->Phi-isrjet->Phi;
			Delta_PT_bl_ISR = bjet1->PT-isrjet->PT;

		      }
		      
		      else {
			//bjet1 is associated to w hadronic
			px_top1=bjet1->PT*cos(bjet1->Phi)+(wjet2->PT*cos(wjet2->Phi))+(wjet1->PT*cos(wjet1->Phi));
			py_top1=bjet1->PT*sin(bjet1->Phi)+(wjet2->PT*sin(wjet2->Phi))+(wjet1->PT*sin(wjet1->Phi));
			
			px_top2=bjet2->PT*cos(bjet2->Phi)+px_lepton;
			py_top2 =bjet2->PT*sin(bjet2->Phi)+py_lepton;
			
			pb1[1]=bjet2->PT*cos(bjet2->Phi);			      
			pb1[2]=bjet2->PT*sin(bjet2->Phi);				
			pb1[3]=bjet2->PT/sin(theta_bjet2);
			pb1[0]=sqrt(bjet2->PT*bjet2->PT+pb1[3]*pb1[3]+bjet2->Mass);
			
			pb2[1]=bjet1->PT*cos(bjet1->Phi);			      
			pb2[2]=bjet1->PT*sin(bjet1->Phi);				
			pb2[3]=bjet1->PT/sin(theta_bjet1);
			pb2[0]=sqrt(bjet1->PT*bjet1->PT+pb2[3]*pb2[3]+bjet1->Mass);
			//Rama hadronica
			Delta_eta_bh_ISR = bjet1->Eta-isrjet->Eta;
			Delta_phi_bh_ISR = bjet1->Phi-isrjet->Phi;
			Delta_PT_bh_ISR = bjet1->PT-isrjet->PT;
			//Rama leptonica
			Delta_eta_bl_ISR = bjet2->Eta-isrjet->Eta;
			Delta_phi_bl_ISR = bjet2->Phi-isrjet->Phi;
			Delta_PT_bl_ISR = bjet2->PT-isrjet->PT;
			
		      }
		      
		      h_Delta_eta_bh_ISR->Fill(Delta_eta_bh_ISR);
		      h_Delta_phi_bh_ISR->Fill(Delta_phi_bh_ISR);
		      h_Delta_PT_bh_ISR->Fill(Delta_PT_bh_ISR);
		      h_Delta_eta_bl_ISR->Fill(Delta_eta_bl_ISR);
		      h_Delta_phi_bl_ISR->Fill(Delta_phi_bl_ISR);
		      h_Delta_PT_bl_ISR->Fill(Delta_PT_bl_ISR);
		      
		      double pa[3]={0};
		      double pb[3]={0};
		      double pmiss[3]={0};
		      pa[0] = mtop;
		      pa[1] = px_top1;
		      pa[2] = py_top1;
		      
		      pb[0] = mtop;
		      pb[1] = px_top2;
		      pb[2] = py_top2;
		      
		      pmiss[0]=0;
		      pmiss[1]= METpointer->MET*cos(METpointer->Phi);
		      pmiss[2]= METpointer->MET*sin(METpointer->Phi);
		      
		      // Calcula MT2
		      mt2_bisect::mt2 mt2_event;
		      mt2_event.set_momenta(pa,pb,pmiss);
		      mt2_event.set_mn(mn);
		      double value_mt2=mt2_event.get_mt2();
		      histMT2->Fill(value_mt2);
		      
		      // Calcula MT2W
		      mt2w_bisect::mt2w mt2w_event;
		      mt2w_event.set_momenta(pl,pb1,pb2,pmiss);
		      double value_mt2w=mt2w_event.get_mt2w();
		      histMT2W->Fill(value_mt2w);
		      
		      //Calcula MT2BL
		      mt2bl_bisect::mt2bl mt2bl_event;
		      mt2bl_event.set_momenta(pl,pb1,pb2,pmiss);
		      double value_mt2bl=mt2bl_event.get_mt2bl();
		      histMT2BL->Fill(value_mt2bl);

		      //Calculo de MTB square

		      double MTB_Square;
		      double theta_bl;
		      double theta_MissingET;
		      double MT_square_MTB_Square;
		      
		      theta_MissingET  = 2*TMath::ATan(TMath::Exp(-METpointer->Eta));
		      theta_bl = TMath::ACos(pz_lepton+pb2z/TMath::Power(pb2x+px_lepton,2)+TMath::Power(pb2y+py_lepton,2)+TMath::Power(pb2z+pz_lepton,2));
		      MTB_Square = 2*METpointer->MET*(eb2+E_lepton-TMath::Sqrt(TMath::Power(pb2x+px_lepton,2)+TMath::Power(pb2y+py_lepton,2)+TMath::Power(pb2z+pz_lepton,2))*TMath::Cos(theta_MissingET-theta_bl));
		      
		      histMTB_Square->Fill(MTB_Square);
		      MT_square_MTB_Square = transverse_mass_square_Muon*MTB_Square;
		      hist_MT_square_MTB_square->Fill(MT_square_MTB_Square);
		      
		      //Subrutina de análisis para los cortes en el PT del ISR y del lepton
		      
		      if (((esolo && Electronpointer->PT>20 && transverse_mass_Electron>80.0&&!(musolo))||
			   (musolo && Muonpointer->PT>20 && transverse_mass_Muon>80.0&&!(esolo)))&&(value_mt2>0.0)) {
			cut_cntr[7]++;//Numero de eventos despues del corte en el lepton en PT y MT 
			rm = METpointer->MET/isrjet->PT;
			rmphi = TMath::Sqrt((rm-0.5)*(rm-0.5)+(Delta_Phi_ISR_MissingET-pi)*(Delta_Phi_ISR_MissingET-pi));
			h_RM_Delta_Phi_ISR_MissingET->Fill(rm,Delta_Phi_ISR_MissingET);
			h_RM->Fill(rm);
			h_RM_PHI->Fill(rmphi);
			
			if((Delta_Phi_ISR_MissingET>0.0)&&(Delta_Phi_ISR_MissingET<6.3)){
			  
			  cut_cntr[8]++; //Numero de eventos despues del corte en el DeltaPhi del ISR y del MissingET
			  if (isrjet->PT>200){
			    for (int jj=0;jj<30;jj++){
			      if(METpointer->MET>10*jj){MET_cntr[jj]++;}
			    }
			    int ij=0;
			    double rmphi_cut=2.5;
			    while(rmphi_cut>0.1){
			      if(rmphi<rmphi_cut)rmphi_cntr[ij]++;
			      rmphi_cut-=0.1;
			      ij++;
			    }
			    ij=0;
			    double rmcut=2.5;
			    while(rmcut>0.1){
			      if(rm<rmcut)rm_cntr[ij]++;
			      rmcut-=0.1;
			      ij++;
			    }
			    ij=0;
			    double mt2cut=350;
			    while(mt2cut>180){
			      if(value_mt2<mt2cut)mt2_cntr[ij]++;
			      mt2cut-=10.0;
			      ij++;
			    }
			    
			    

			  }
			  
			  if(METpointer->MET>100){
			    cut_cntr[9]++; //Numero de eventos despues del corte en el MissingET > 100 GeV
			    
			    if (isrjet->PT>100){
				  cut_cntr[10]++;//Numero de eventos despues del corte en el PT del ISR > 100 GeV
				  
				  h_Delta_eta_bh_ISR_corte_100_GeV->Fill(Delta_eta_bh_ISR);
				  h_Delta_phi_bh_ISR_corte_100_GeV->Fill(Delta_phi_bh_ISR);
				  h_Delta_PT_bh_ISR_corte_100_GeV->Fill(Delta_PT_bh_ISR);
				  h_Delta_eta_bl_ISR_corte_100_GeV->Fill(Delta_eta_bl_ISR);
				  h_Delta_phi_bl_ISR_corte_100_GeV->Fill(Delta_phi_bl_ISR);
				  h_Delta_PT_bl_ISR_corte_100_GeV->Fill(Delta_PT_bl_ISR);
				  h_MET_corte_100_GeV->Fill(MET);
				  h_MT_Muon_corte_100_GeV->Fill(transverse_mass_Muon);
				  h_MT_Electron_corte_100_GeV->Fill(transverse_mass_Electron);
				  h_MT_Muon_square_corte_100_GeV->Fill(transverse_mass_square_Muon);
				  h_MT_Electron_square_corte_100_GeV->Fill(transverse_mass_square_Electron);
				  h_Cociente_PT_Muon_PT_ISR_corte_100_GeV->Fill(Cociente_PT_Muon_PT_ISR);
				  h_Cociente_PT_Electron_PT_ISR_corte_100_GeV->Fill(Cociente_PT_Electron_PT_ISR);
				  h_Delta_Eta_Muon_ISR_corte_100_GeV->Fill(Delta_Eta_Muon_ISR);
				  h_Delta_Phi_Muon_ISR_corte_100_GeV->Fill(Delta_Phi_Muon_ISR);
				  h_Delta_Eta_Electron_ISR_corte_100_GeV->Fill(Delta_Eta_Electron_ISR);
				  h_Delta_Phi_Electron_ISR_corte_100_GeV->Fill(Delta_Phi_Electron_ISR);
				  h_Cociente_MissingET_PT_ISR_corte_100_GeV->Fill(Cociente_MissingET_PT_ISR);
				  h_Delta_Eta_ISR_MissingET_corte_100_GeV->Fill(Delta_Eta_ISR_MissingET);
				  h_Delta_Phi_ISR_MissingET_corte_100_GeV->Fill(Delta_Phi_ISR_MissingET);
				  h_PT_Muon_MissingET_corte_100_GeV->Fill(PT_Muon,MissingET);
				  h_PT_Electron_MissingET_corte_100_GeV->Fill(PT_Electron,MissingET);
				  h_PT_Muon_Delta_Eta_Muon_ISR_corte_100_GeV->Fill(PT_Muon,Delta_Eta_Muon_ISR);
				  h_PT_Electron_Delta_Eta_Electron_ISR_100_GeV->Fill(PT_Electron,Delta_Eta_Electron_ISR);
				  h_PT_Muon_PT_ISR_corte_100_GeV->Fill(PT_Muon,PT_ISR);
				  h_PT_Electron_PT_ISR_corte_100_GeV->Fill(PT_Electron,PT_ISR);
				  histMT2_corte_100_GeV->Fill(value_mt2);	
				  histMT2W_corte_100_GeV->Fill(value_mt2w);
				  histMT2BL_corte_100_GeV->Fill(value_mt2bl);
				  histMTB_Square_corte_100_GeV->Fill(MTB_Square);
				  h_RM_Delta_Phi_ISR_MissingET_corte_100_GeV->Fill(rm,Delta_Phi_ISR_MissingET);
				  h_RM_corte_100_GeV->Fill(rm);
				  h_RM_PHI_corte_100_GeV->Fill(rmphi);	
				  
				}
				if (isrjet->PT>300) {
				  cut_cntr[11]++; //Numero de eventos despues del corte en el PT del ISR > 300 GeV
				  h_Delta_eta_bh_ISR_corte_300_GeV->Fill(Delta_eta_bh_ISR);
				  h_Delta_phi_bh_ISR_corte_300_GeV->Fill(Delta_phi_bh_ISR);
				  h_Delta_PT_bh_ISR_corte_300_GeV->Fill(Delta_PT_bh_ISR);
				  h_Delta_eta_bl_ISR_corte_300_GeV->Fill(Delta_eta_bl_ISR);
				  h_Delta_phi_bl_ISR_corte_300_GeV->Fill(Delta_phi_bl_ISR);
				  h_Delta_PT_bl_ISR_corte_300_GeV->Fill(Delta_PT_bl_ISR);
				  h_MET_corte_300_GeV->Fill(MET);
				  h_MT_Muon_corte_300_GeV->Fill(transverse_mass_Muon);
				  h_MT_Electron_corte_300_GeV->Fill(transverse_mass_Electron);
				  h_MT_Muon_square_corte_300_GeV->Fill(transverse_mass_square_Muon);
				  h_MT_Electron_square_corte_300_GeV->Fill(transverse_mass_square_Electron);
				  h_Cociente_PT_Muon_PT_ISR_corte_300_GeV->Fill(Cociente_PT_Muon_PT_ISR);
				  h_Cociente_PT_Electron_PT_ISR_corte_300_GeV->Fill(Cociente_PT_Electron_PT_ISR);
				  h_Delta_Eta_Muon_ISR_corte_300_GeV->Fill(Delta_Eta_Muon_ISR);
				  h_Delta_Phi_Muon_ISR_corte_300_GeV->Fill(Delta_Phi_Muon_ISR);
				  h_Delta_Eta_Electron_ISR_corte_300_GeV->Fill(Delta_Eta_Electron_ISR);
				  h_Delta_Phi_Electron_ISR_corte_300_GeV->Fill(Delta_Phi_Electron_ISR);
				  h_Cociente_MissingET_PT_ISR_corte_300_GeV->Fill(Cociente_MissingET_PT_ISR);
				  h_Delta_Eta_ISR_MissingET_corte_300_GeV->Fill(Delta_Eta_ISR_MissingET);
				  h_Delta_Phi_ISR_MissingET_corte_300_GeV->Fill(Delta_Phi_ISR_MissingET);
				  h_PT_Muon_MissingET_corte_300_GeV->Fill(PT_Muon,MissingET);
				  h_PT_Electron_MissingET_corte_300_GeV->Fill(PT_Electron,MissingET);
				  h_PT_Muon_Delta_Eta_Muon_ISR_corte_300_GeV->Fill(PT_Muon,Delta_Eta_Muon_ISR);
				  h_PT_Electron_Delta_Eta_Electron_ISR_300_GeV->Fill(PT_Electron,Delta_Eta_Electron_ISR);
				  h_PT_Muon_PT_ISR_corte_300_GeV->Fill(PT_Muon,PT_ISR);
				  h_PT_Electron_PT_ISR_corte_300_GeV->Fill(PT_Electron,PT_ISR);
				  histMT2_corte_300_GeV->Fill(value_mt2);
				  histMT2W_corte_300_GeV->Fill(value_mt2w);
				  histMT2BL_corte_300_GeV->Fill(value_mt2bl);
				  histMTB_Square_corte_300_GeV->Fill(MTB_Square);
				  h_RM_Delta_Phi_ISR_MissingET_corte_300_GeV->Fill(rm,Delta_Phi_ISR_MissingET);
				  h_RM_corte_300_GeV->Fill(rm);
				  h_RM_PHI_corte_300_GeV->Fill(rmphi);	
					
				  
				}
				if (isrjet->PT>500){
				  cut_cntr[12]++; //Numero de eventos despues del corte en el PT del ISR > 500 GeV
				  h_Delta_eta_bh_ISR_corte_500_GeV->Fill(Delta_eta_bh_ISR);
				  h_Delta_phi_bh_ISR_corte_500_GeV->Fill(Delta_phi_bh_ISR);
				  h_Delta_PT_bh_ISR_corte_500_GeV->Fill(Delta_PT_bh_ISR);
				  h_Delta_eta_bl_ISR_corte_500_GeV->Fill(Delta_eta_bl_ISR);
				  h_Delta_phi_bl_ISR_corte_500_GeV->Fill(Delta_phi_bl_ISR);
				  h_Delta_PT_bl_ISR_corte_500_GeV->Fill(Delta_PT_bl_ISR);
				  h_MET_corte_500_GeV->Fill(MET);
				  h_MT_Muon_corte_500_GeV->Fill(transverse_mass_Muon);
				  h_MT_Electron_corte_500_GeV->Fill(transverse_mass_Electron);
				  h_MT_Muon_square_corte_500_GeV->Fill(transverse_mass_square_Muon);
				  h_MT_Electron_square_corte_500_GeV->Fill(transverse_mass_square_Electron);
				  h_Cociente_PT_Muon_PT_ISR_corte_500_GeV->Fill(Cociente_PT_Muon_PT_ISR);
				  h_Cociente_PT_Electron_PT_ISR_corte_500_GeV->Fill(Cociente_PT_Electron_PT_ISR);
				  h_Delta_Eta_Muon_ISR_corte_500_GeV->Fill(Delta_Eta_Muon_ISR);
				  h_Delta_Phi_Muon_ISR_corte_500_GeV->Fill(Delta_Phi_Muon_ISR);
				  h_Delta_Eta_Electron_ISR_corte_500_GeV->Fill(Delta_Eta_Electron_ISR);
				  h_Delta_Phi_Electron_ISR_corte_500_GeV->Fill(Delta_Phi_Electron_ISR);
				  h_Cociente_MissingET_PT_ISR_corte_500_GeV->Fill(Cociente_MissingET_PT_ISR);
				  h_Delta_Eta_ISR_MissingET_corte_500_GeV->Fill(Delta_Eta_ISR_MissingET);
				  h_Delta_Phi_ISR_MissingET_corte_500_GeV->Fill(Delta_Phi_ISR_MissingET);
				  h_PT_Muon_MissingET_corte_500_GeV->Fill(PT_Muon,MissingET);
				  h_PT_Electron_MissingET_corte_500_GeV->Fill(PT_Electron,MissingET);
				  h_PT_Muon_Delta_Eta_Muon_ISR_corte_500_GeV->Fill(PT_Muon,Delta_Eta_Muon_ISR);
				  h_PT_Electron_Delta_Eta_Electron_ISR_500_GeV->Fill(PT_Electron,Delta_Eta_Electron_ISR);
				  h_PT_Muon_PT_ISR_corte_500_GeV->Fill(PT_Muon,PT_ISR);
				  h_PT_Electron_PT_ISR_corte_500_GeV->Fill(PT_Electron,PT_ISR);
				  histMT2_corte_500_GeV->Fill(value_mt2);
				  histMT2W_corte_500_GeV->Fill(value_mt2w);
				  histMT2BL_corte_500_GeV->Fill(value_mt2bl);	
				  histMTB_Square_corte_500_GeV->Fill(MTB_Square);	
				  h_RM_Delta_Phi_ISR_MissingET_corte_500_GeV->Fill(rm,Delta_Phi_ISR_MissingET);
				  h_RM_corte_500_GeV->Fill(rm);
				  h_RM_PHI_corte_500_GeV->Fill(rmphi);	
				
	
				}
				}
			      }
		      }
		      
		      }
		      
		      }
		    }
		  }	      //loop ISR_found
		  
		  
		  // Jet with greatest PT
		  if (posLeadingPT != -1){
		    h_leading_PT->Fill(PT_max);
		    if(posLeadingPT == ISR_jets[entry]) ISR_greatest_PT++;

		    currentJet = (Jet*) branchJet->At(posLeadingPT);
		    vect_currentJet->SetPtEtaPhiM(currentJet->PT,currentJet->Eta,currentJet->Phi,currentJet->Mass);
		    delta_phi = deltaAng(vect_currentJet->Phi(), METpointer->Phi);
		    MT_leading_jet = sqrt(2*vect_currentJet->Pt()*MET*(1-cos(delta_phi)));
		    h_leading_MT->Fill(MT_leading_jet);
			  
		    h_leading_Eta->Fill(vect_currentJet->Eta());
		    h_leading_DPhi_MET->Fill(delta_phi);
		    
		    h2_leading_PTEta->Fill(vect_currentJet->Pt(),vect_currentJet->Eta());
		  }
		  
		  // HT
		  if (1 < HT_R1/HT || 1 < HT_R2/HT){
		    cout << "Error en el evento: " << entry << endl;
		    cout << "HT: " << HT << "\tHT_R1: " << HT_R1 << "\tHT_R2: " << HT_R2 << endl;
		    return 1;
		  }
		  
		  h_HT->Fill(HT);
		  h_HT_R1->Fill(HT_R1/HT);
		  h_HT_R2->Fill(HT_R2/HT);
		  h_BTags_per_Event->Fill(number_Btags);
		  
		} //end loop over events
		
		cout<<"progress = 100%\t";
		cout<< "Time :"<< (clock()-initialTime)/double_t(CLOCKS_PER_SEC)<<"s"<<endl;
		cout<< "Percentage of events where the ISR jet is the jet with greatest PT: " << (Double_t) (ISR_greatest_PT*100)/numberOfEntries << "%\n";
		cout<< "Percentage of events where the ISR jet is tagged as Bjet: " << (Double_t) (ISR_Btags*100)/numberOfEntries << "%\n";
		
		
	} // End run's for cicle
	cout<<" ****************************************"<<endl;
	cout<<" # events with 5 jets ="<<n5jets<<endl;
	cout<<" ****************************************"<<endl;
	
	cout<<" **********EVENT COUNTER******"<<endl;
	cout<<" ****************************************"<<endl;
	
	cout<<"Total events processed                                  ="<<cut_cntr[0]<<endl;
	cout<<"No. events with muon or electron                        ="<<cut_cntr[1]<<endl;
	cout<<"No. events with ISR + isolated lepton                   ="<<cut_cntr[2]<<endl;
	cout<<"No. events with 1 bjet                                  ="<<cut_cntr[3]<<endl;
	cout<<"No. events with 2 bjet                                  ="<<cut_cntr[4]<<endl;
	cout<<"No. events with 2 wjet                                  ="<<cut_cntr[5]<<endl;
	cout<<"No. events with 2 bjet and cut PT                       ="<<cut_cntr[6]<<endl;
	cout<<"No. events after lepton PT and MT cut                   ="<<cut_cntr[7]<<endl;
	cout<<"No. events after DeltaPhi cut between ISR and MissingET ="<<cut_cntr[8]<<endl;
	cout<<"No. events after MissingET cut > 100 GeV                ="<<cut_cntr[9]<<endl;
	cout<<"No. events after ISR_PT>100 GeV                         ="<<cut_cntr[10]<<endl;
	cout<<"No. events after ISR_PT>300 GeV                         ="<<cut_cntr[11]<<endl;
	cout<<"No. events after ISR_PT>500 GeV                         ="<<cut_cntr[12]<<endl;

 	cout<<" ****************************************"<<endl;
	cout<<" "<<endl;
	cout<<" "<<endl;
	cout<<" ****************************************"<<endl;
	cout<<" *************MET SCAN************"<<endl;
	cout<<" ****************************************"<<endl;
	for (int jj=0;jj<30;jj++){
	  cout<<"# events with MET > "<<jj*10<<" GeV  = "<<MET_cntr[jj]<<endl;
	}
	cout<<" ****************************************"<<endl;
	cout<<" *************RMPHI SCAN************"<<endl;
	cout<<" ****************************************"<<endl;
	for (int jj=0;jj<24;jj++){
	  cout<<"# events with RMPHI<  "<<(25-jj)*0.1<< "   = "<<rmphi_cntr[jj]<<endl;
	}
	cout<<" ****************************************"<<endl;
	cout<<" *************RM SCAN************"<<endl;
	cout<<" ****************************************"<<endl;
	for (int jj=0;jj<24;jj++){
	  cout<<"# events with RM<  "<<(25-jj)*0.1<< "   = "<<rm_cntr[jj]<<endl;
	}

	cout<<" ****************************************"<<endl;
	cout<<" *************MT2 SCAN************"<<endl;
	cout<<" ****************************************"<<endl;
	for (int jj=0;jj<17;jj++){
	  cout<<"# events with MT2>  "<<(18+jj)*10.0<< "   = "<<mt2_cntr[jj]<<endl;
	}


	cout<<" ****************************************"<<endl;

	string hfile_name_str = head_folder_results + "histos_Stops_200_N1_20.root";
	//string hfile_name_str = head_folder_results + "histos_Stops_200_N1_35.root";
	//string hfile_name_str = head_folder_results + "histos_Stops_300_N1_135.root";
	//string hfile_name_str = head_folder_results + "histos_Stops_300_N1_120.root";
	//string hfile_name_str = head_folder_results + "histos_Stops_400_N1_220.root";
	//string hfile_name_str = head_folder_results + "histos_Stops_400_N1_235.root";
	//string hfile_name_str = head_folder_results + "histos_Stops_500_N1_320.root";
	//string hfile_name_str = head_folder_results + "histos.root";
	//string hfile_name_str = head_folder_results + "histos_Top_semileptónico.root";	
	//string hfile_name_str = head_folder_results + "histos_Top_leptónico.root";


	Char_t *hfile_name = (Char_t *) hfile_name_str.c_str();

	TFile* hfile = new TFile(hfile_name, "RECREATE");
	h_jet_DPhi_MET->Write();
	h_jet_Eta->Write();
	h_MuonPT->Write();
	h_jet_PT->Write();
	h_jet_Phi->Write();
	h_jet_MT->Write();
	h_MT_Muon->Write();
	h_MT_Muon_corte_100_GeV->Write();
	h_MT_Muon_corte_300_GeV->Write();
	h_MT_Muon_corte_500_GeV->Write();
	h_MT_Electron->Write();
	h_MT_Electron_corte_100_GeV->Write();
	h_MT_Electron_corte_300_GeV->Write();
	h_MT_Electron_corte_500_GeV->Write();
	h_MT_square_Muon->Write();
	h_MT_Muon_square_corte_100_GeV->Write();
	h_MT_Muon_square_corte_300_GeV->Write();
	h_MT_Muon_square_corte_500_GeV->Write();
	h_MT_square_Electron->Write();
	h_MT_Electron_square_corte_100_GeV->Write();
	h_MT_Electron_square_corte_300_GeV->Write();
	h_MT_Electron_square_corte_500_GeV->Write();
	h_ISR_MET_HT->Write();
	h_ISR_MET_HT1->Write();
	h_ISR_MET_HT2->Write();
	h_ISR_MET_HT3->Write();
	h_jet_Delta_PT->Write();
	h_jet_Delta_Eta->Write();
	h_jet_DPhi_MET_other->Write();
	h_jet_PT_HT->Write();
	h_jet_multiplicity->Write();
	h_jet_PT_over_PT_others->Write();
	h_jet_Eta_over_Eta_others->Write();
	h_jet_DeltaR->Write();
	h_jet_DPhi_over_Phi_others->Write();
	h_jet_Delta_Eta_leading->Write();
	h_jet_Delta_PT_leading->Write();
	//Se escriben los histogramas de MT2, MT2W, MT2BL, MTB Square
	histMT2->Write();
	histMT2_corte_100_GeV->Write();
	histMT2_corte_300_GeV->Write();
	histMT2_corte_500_GeV->Write();
	histMT2W->Write();
	histMT2W_corte_100_GeV->Write();
	histMT2W_corte_300_GeV->Write();
	histMT2W_corte_500_GeV->Write();
	histMT2BL->Write();
	histMT2BL_corte_100_GeV->Write();
	histMT2BL_corte_300_GeV->Write();
	histMT2BL_corte_500_GeV->Write();
	histMTB_Square->Write();
	histMTB_Square_corte_100_GeV->Write();
	histMTB_Square_corte_300_GeV->Write();
	histMTB_Square_corte_500_GeV->Write();
	hist_MT_square_MTB_square->Write();
	//Correlación entre variables
	h_Delta_eta_bh_ISR->Write();
	h_Delta_phi_bh_ISR->Write();
	h_Delta_PT_bh_ISR->Write();
	h_Delta_eta_bl_ISR->Write();
	h_Delta_phi_bl_ISR->Write();
	h_Delta_PT_bl_ISR->Write();
	h_Delta_eta_bh_ISR_corte_100_GeV->Write();
	h_Delta_phi_bh_ISR_corte_100_GeV->Write();
	h_Delta_PT_bh_ISR_corte_100_GeV->Write();
	h_Delta_eta_bl_ISR_corte_100_GeV->Write();
	h_Delta_phi_bl_ISR_corte_100_GeV->Write();
	h_Delta_PT_bl_ISR_corte_100_GeV->Write();
	h_Delta_eta_bh_ISR_corte_300_GeV->Write();
	h_Delta_phi_bh_ISR_corte_300_GeV->Write();
	h_Delta_PT_bh_ISR_corte_300_GeV->Write();
	h_Delta_eta_bl_ISR_corte_300_GeV->Write();
	h_Delta_phi_bl_ISR_corte_300_GeV->Write();
	h_Delta_PT_bl_ISR_corte_300_GeV->Write();
	h_Delta_eta_bh_ISR_corte_500_GeV->Write();
	h_Delta_phi_bh_ISR_corte_500_GeV->Write();
	h_Delta_PT_bh_ISR_corte_500_GeV->Write();
	h_Delta_eta_bl_ISR_corte_500_GeV->Write();
	h_Delta_phi_bl_ISR_corte_500_GeV->Write();
	h_Delta_PT_bl_ISR_corte_500_GeV->Write();
	h_Cociente_PT_Muon_PT_ISR->Write();
	h_Cociente_PT_Muon_PT_ISR_corte_100_GeV->Write();
	h_Cociente_PT_Muon_PT_ISR_corte_300_GeV->Write();
	h_Cociente_PT_Muon_PT_ISR_corte_500_GeV->Write();
	h_Cociente_PT_Electron_PT_ISR->Write();
	h_Cociente_PT_Electron_PT_ISR_corte_100_GeV->Write();
	h_Cociente_PT_Electron_PT_ISR_corte_300_GeV->Write();
	h_Cociente_PT_Electron_PT_ISR_corte_500_GeV->Write();
	h_Delta_Eta_Muon_ISR->Write();
	h_Delta_Eta_Muon_ISR_corte_100_GeV->Write();
	h_Delta_Eta_Muon_ISR_corte_300_GeV->Write();
	h_Delta_Eta_Muon_ISR_corte_500_GeV->Write();
	h_Delta_Phi_Muon_ISR->Write();
	h_Delta_Phi_Muon_ISR_corte_100_GeV->Write();
	h_Delta_Phi_Muon_ISR_corte_300_GeV->Write();
	h_Delta_Phi_Muon_ISR_corte_500_GeV->Write();
	h_Delta_Eta_Electron_ISR->Write();
	h_Delta_Eta_Electron_ISR_corte_100_GeV->Write();
	h_Delta_Eta_Electron_ISR_corte_300_GeV->Write();
	h_Delta_Eta_Electron_ISR_corte_500_GeV->Write();
	h_Delta_Phi_Electron_ISR->Write();
	h_Delta_Phi_Electron_ISR_corte_100_GeV->Write();
	h_Delta_Phi_Electron_ISR_corte_300_GeV->Write();
	h_Delta_Phi_Electron_ISR_corte_500_GeV->Write();
	h_Cociente_MissingET_PT_ISR->Write();
	h_Cociente_MissingET_PT_ISR_corte_100_GeV->Write();
	h_Cociente_MissingET_PT_ISR_corte_300_GeV->Write();
	h_Cociente_MissingET_PT_ISR_corte_500_GeV->Write();
	h_Delta_Eta_ISR_MissingET->Write();
	h_Delta_Eta_ISR_MissingET_corte_100_GeV->Write();
	h_Delta_Eta_ISR_MissingET_corte_300_GeV->Write();
	h_Delta_Eta_ISR_MissingET_corte_500_GeV->Write();
	h_Delta_Phi_ISR_MissingET->Write();
	h_Delta_Phi_ISR_MissingET_corte_100_GeV->Write();
	h_Delta_Phi_ISR_MissingET_corte_300_GeV->Write();
	h_Delta_Phi_ISR_MissingET_corte_500_GeV->Write();
	h_PT_Muon_MissingET->Write();
	h_PT_Muon_MissingET_corte_100_GeV->Write();
	h_PT_Muon_MissingET_corte_300_GeV->Write();
	h_PT_Muon_MissingET_corte_500_GeV->Write();
	h_PT_Electron_MissingET->Write();
	h_PT_Electron_MissingET_corte_100_GeV->Write();
	h_PT_Electron_MissingET_corte_300_GeV->Write();
	h_PT_Electron_MissingET_corte_500_GeV->Write();
	h_PT_Muon_Delta_Eta_Muon_ISR->Write();
	h_PT_Muon_Delta_Eta_Muon_ISR_corte_100_GeV->Write();
	h_PT_Muon_Delta_Eta_Muon_ISR_corte_300_GeV->Write();
	h_PT_Muon_Delta_Eta_Muon_ISR_corte_500_GeV->Write();
	h_PT_Electron_Delta_Eta_Electron_ISR->Write();
	h_PT_Electron_Delta_Eta_Electron_ISR_100_GeV->Write();
	h_PT_Electron_Delta_Eta_Electron_ISR_300_GeV->Write();
	h_PT_Electron_Delta_Eta_Electron_ISR_500_GeV->Write();
	h_Delta_Eta_Muon_ISR_Delta_Eta_ISR_MissingET->Write();
	h_Delta_Eta_Muon_ISR_Delta_Eta_ISR_MissingET_corte_100_GeV->Write();
	h_Delta_Eta_Muon_ISR_Delta_Eta_ISR_MissingET_corte_300_GeV->Write();
	h_Delta_Eta_Muon_ISR_Delta_Eta_ISR_MissingET_corte_500_GeV->Write();
	h_Delta_Eta_Electron_ISR_Delta_Eta_ISR_MissingET->Write();
	h_Delta_Eta_Electron_ISR_Delta_Eta_ISR_MissingET_corte_100_GeV->Write();
	h_Delta_Eta_Electron_ISR_Delta_Eta_ISR_MissingET_corte_300_GeV->Write();
	h_Delta_Eta_Electron_ISR_Delta_Eta_ISR_MissingET_corte_500_GeV->Write();
	h_PT_Muon_PT_ISR->Write();
	h_PT_Muon_PT_ISR_corte_100_GeV->Write();
	h_PT_Muon_PT_ISR_corte_300_GeV->Write();
	h_PT_Muon_PT_ISR_corte_500_GeV->Write();
	h_PT_Electron_PT_ISR->Write();
	h_PT_Electron_PT_ISR_corte_100_GeV->Write();
	h_PT_Electron_PT_ISR_corte_300_GeV->Write();
	h_PT_Electron_PT_ISR_corte_500_GeV->Write();
	h_RM_Delta_Phi_ISR_MissingET->Write();
	h_RM->Write();
	h_RM_PHI->Write();
	h_RM_Delta_Phi_ISR_MissingET_corte_100_GeV->Write();
	h_RM_corte_100_GeV->Write();
	h_RM_PHI_corte_100_GeV->Write();
	h_RM_Delta_Phi_ISR_MissingET_corte_300_GeV->Write();
	h_RM_corte_300_GeV->Write();
	h_RM_PHI_corte_300_GeV->Write();
	h_RM_Delta_Phi_ISR_MissingET_corte_500_GeV->Write();
	h_RM_corte_500_GeV->Write();
	h_RM_PHI_corte_500_GeV->Write();

	//Final de correlaciones	
	h_ISR_DPhi_MET->Write();
	h_ISR_Eta->Write();
	h_ISR_PT->Write();
	h_ISR_Phi->Write();
	//h_ISR_MT->Write();
	h_ISR_Delta_PT->Write();
	h_ISR_Delta_Eta->Write();
	h_ISR_DPhi_MET_other->Write();
	h_ISR_PT_HT->Write();
	h_ISR_multiplicity->Write();
	h_ISR_PT_over_PT_others->Write();
	h_ISR_Eta_over_Eta_others->Write();
	h_ISR_DeltaR->Write();
	h_ISR_DPhi_over_Phi_others->Write();
	h_ISR_Delta_Eta_leading->Write();
	h_ISR_Delta_PT_leading->Write();
	
	h_MET->Write();
	h_MET_corte_100_GeV->Write();
	h_MET_corte_300_GeV->Write();
	h_MET_corte_500_GeV->Write();
	
	h_leading_MT->Write();
	h_leading_PT->Write();
	h_leading_Eta->Write();
	h_leading_DPhi_MET->Write();

	h_HT->Write();
	h_HT_R1->Write();
	h_HT_R2->Write();

	h_numberJet->Write();
	
	h_BTag->Write();
	h_BTag_PT->Write();
	h_BTag_Eta->Write();
	h_BTag_DPhi_MET->Write();
	h_BTags_per_Event->Write();
	
	h2_ISR_PTEta->Write();
	h2_jet_PTEta->Write();
	h2_dif_PTEta->Add(h2_ISR_PTEta,h2_jet_PTEta,1,-1);
	h2_dif_PTEta->Write();

	h2_dif_lead_PTEta->Add(h2_ISR_PTEta,h2_leading_PTEta,1,-1);
	h2_dif_lead_PTEta->Write();
	
	{
	  string salida_str = head_folder_results + "Eta";
	  Char_t *salida = (Char_t *) salida_str.c_str();
	  TCanvas *C = new TCanvas(salida,"Pseudorapidity",1280,720);
	  Present(h_ISR_Eta,h_jet_Eta,C,1,"h","Num. Jets / Total",122);
	  C->Write();
	  C->Close();
	  
	  salida_str = head_folder_results + "Eta ISR vs BTag";
	  salida = (Char_t *) salida_str.c_str();
	  C = new TCanvas(salida,"Pseudorapidity ISR vs BTag",1280,720);
	  Present(h_ISR_Eta,h_BTag_Eta,C,1,"h","Num. Jets / Total",122);
	  C->Write();
	  C->Close();
		
	  salida_str = head_folder_results + "Eta ISR vs Leading";
	  salida = (Char_t *) salida_str.c_str();
	  C = new TCanvas(salida,"Pseudorapidity ISR vs Leading",1280,720);
	  Present(h_ISR_Eta,h_leading_Eta,C,1,"h","Num. Jets / Total",122);
	  C->Write();
	  C->Close();
		
	  salida_str = head_folder_results + "Transverse momentum";
	  salida = (Char_t *) salida_str.c_str();
	  C = new TCanvas(salida,"Transverse momentum",1280,720);
	  Present(h_ISR_PT,h_jet_PT,C,2,"PT [GeV]","Num. Jets / Total");
	  C->Write();
	  C->Close();
	  
	  salida_str = head_folder_results + "Transverse momentum ISR vs Leading";
	  salida = (Char_t *) salida_str.c_str();
	  C = new TCanvas(salida,"Transverse momentum ISR vs Leading",1280,720);
	  Present(h_ISR_PT,h_leading_PT,C,2,"PT [GeV]","Num. Jets / Total");
	  C->Write();
	  C->Close();
		
	  salida_str = head_folder_results + "Transverse momentum ISR vs B_Tag";
	  salida = (Char_t *) salida_str.c_str();
	  C = new TCanvas(salida,"Transverse momentum ISR vs B_Tag",1280,720);
	  Present(h_ISR_PT,h_BTag_PT,C,2,"PT [GeV]","Num. Jets / Total");
	  C->Write();
	  C->Close();
	  
	  salida_str = head_folder_results + "Transverse momentum ISR, B_Tag, Leading";
	  salida = (Char_t *) salida_str.c_str();
	  C = new TCanvas(salida,"Transverse momentum ISR, B_Tag, Leading",1280,720);
	  Present_3(h_ISR_PT,h_BTag_PT,h_leading_PT,C,2,"PT [GeV]","Num. Jets / Total");
	  C->Write();
	  C->Close();
	  
	  salida_str = head_folder_results + "Transverse momentum ISR, B_Tag, Leading LOG";
	  salida = (Char_t *) salida_str.c_str();
	  C = new TCanvas(salida,"Transverse momentum ISR, B_Tag, Leading LOG",1280,720);
	  Present_3(h_ISR_PT,h_BTag_PT,h_leading_PT,C,2,"PT [GeV]","Num. Jets / Total",12,12,true);
	  C->Write();
	  C->Close();
	  
	  salida_str = head_folder_results + "Transverse mass Leading vs ISR Jet";
	  salida = (Char_t *) salida_str.c_str();
	  C = new TCanvas(salida,"Transverse mass Leading vs ISR Jet",1280,720);
	  Present(h_MT_Muon,h_leading_MT,C,2,"MT [GeV]","Num. Jets / Total");
	  C->Write();
	  C->Close();
		
	  salida_str = head_folder_results + "Transverse mass ISR vs Jet";
	  salida = (Char_t *) salida_str.c_str();
	  C = new TCanvas(salida,"Transverse mass ISR vs Jet",1280,720);
	  Present(h_MT_Muon,h_jet_MT,C,2,"MT [GeV]","Num. Jets / Total");
	  C->Write();
	  C->Close();

	  salida_str = head_folder_results + "Phi";
	  salida = (Char_t *) salida_str.c_str();
	  C = new TCanvas(salida,"Phi",1280,720);
	  Present(h_ISR_Phi,h_jet_Phi,C,3,"f","Num. Jets / Total",122);
	  C->Write();
	  C->Close();
		
	  salida_str = head_folder_results + "Delta Phi - Jet - MET";
	  salida = (Char_t *) salida_str.c_str();
	  C = new TCanvas(salida,"Delta Phi - Jet - MET",1280,720);
	  Present(h_ISR_DPhi_MET,h_jet_DPhi_MET,C,3,"Df","Num. Jets / Total",122);
	  C->Write();
	  C->Close();

	  salida_str = head_folder_results + "Delta Phi - Jet - MET - Btag";
	  salida = (Char_t *) salida_str.c_str();
	  C = new TCanvas(salida,"Delta Phi - Jet - MET - Btag",1280,720);
	  Present(h_ISR_DPhi_MET,h_BTag_DPhi_MET,C,3,"Df","Num. Jets / Total",122);
	  C->Write();
	  C->Close();

	  salida_str = head_folder_results + "Delta Phi - Jet - MET - leading";
	  salida = (Char_t *) salida_str.c_str();
	  C = new TCanvas(salida,"Delta Phi - Jet - MET - leading",1280,720);
	  Present(h_ISR_DPhi_MET,h_leading_DPhi_MET,C,1,"Df","Num. Jets / Total",122);
	  C->Write();
	  C->Close();

	  salida_str = head_folder_results + "MET > 120";
	  salida = (Char_t *) salida_str.c_str();
	  C = new TCanvas(salida,"MET > 120",1280,720);
	  Present(h_MET,h_MET_corte_100_GeV,C,2,"MET","Num. Jets / Total");
	  C->Write();
	  C->Close();
		
	  salida_str = head_folder_results + "MET > 200";
	  salida = (Char_t *) salida_str.c_str();
	  C = new TCanvas(salida,"MET > 200",1280,720);
	  Present(h_MET,h_MET_corte_300_GeV,C,2,"MET","Num. Jets / Total");
	  C->Write();
	  C->Close();
	  
	  salida_str = head_folder_results + "MET > 240";
	  salida = (Char_t *) salida_str.c_str();
	  C = new TCanvas(salida,"MET > 240",1280,720);
	  Present(h_MET,h_MET_corte_500_GeV,C,2,"MET","Num. Jets / Total");
	  C->Write();
	  C->Close();
	  
	  salida_str = head_folder_results + "HT ratio comparison";
	  salida = (Char_t *) salida_str.c_str();
	  C = new TCanvas(salida,"HT ratio comparison",1280,720);
	  Present(h_HT_R1,h_HT_R2,C,2,"HT","Num. Jets / Total");
	  C->Write();
	  C->Close();
		
	  salida_str = head_folder_results + "PT vs ETA - ISR";
	  salida = (Char_t *) salida_str.c_str();
	  C = new TCanvas(salida,"PT vs ETA - ISR",1280,720);
	  Plot_Single_2D(h2_ISR_PTEta,C,2, "PT [GeV]", "h", 12, 122);
	  C->Write();
	  C->Close();
		
	  salida_str = head_folder_results + "PT vs ETA - Jet";
	  salida = (Char_t *) salida_str.c_str();
	  C = new TCanvas(salida,"PT vs ETA - Jet",1280,720);
	  Plot_Single_2D(h2_jet_PTEta,C,2, "PT [GeV]", "h", 12, 122);
	  C->Write();
	  C->Close();
		
	  salida_str = head_folder_results + "PT vs ETA - Diff with any jet";
	  salida = (Char_t *) salida_str.c_str();
	  C = new TCanvas(salida,"PT vs ETA - Diff with any jet",1280,720);
	  Plot_Single_2D(h2_dif_PTEta,C,2, "PT [GeV]", "h", 12, 122);
	  C->Write();
	  C->Close();
		
	  salida_str = head_folder_results + "PT vs ETA - leading";
	  salida = (Char_t *) salida_str.c_str();
	  C = new TCanvas(salida,"PT vs ETA - leading",1280,720);
	  Plot_Single_2D(h2_leading_PTEta,C,2, "PT [GeV]", "h", 12, 122);
	  C->Write();
	  C->Close();
		
	  salida_str = head_folder_results + "PT vs ETA - Diff with leading";
	  salida = (Char_t *) salida_str.c_str();
	  C = new TCanvas(salida,"PT vs ETA - Diff with leading",1280,720);
	  Plot_Single_2D(h2_dif_lead_PTEta,C,2, "PT [GeV]", "h", 12, 122);
	  C->Write();
	  C->Close();
	  
	  salida_str = head_folder_results + "HT";
	  salida = (Char_t *) salida_str.c_str();
	  C = new TCanvas(salida,"HT",1280,720);
	  Plot_Single(h_HT,C,2, "HT [GeV]", "Num. Jets / Total", 12, 12);
	  C->Write();
	  C->Close();

	  salida_str = head_folder_results + "Number_of_B_Tags";
	  salida = (Char_t *) salida_str.c_str();
	  C = new TCanvas(salida,"Number of B Tags",1280,720);
	  Plot_Single(h_BTags_per_Event,C,2, "B Tags / event", "Num. Jets / Total", 12, 12);
	  C->Write();
	  C->Close();
	  
	  salida_str = head_folder_results + "Jet_multiplitcity";
	  salida = (Char_t *) salida_str.c_str();
	  C = new TCanvas(salida,"Jet multiplicity",1280,720);
	  Present(h_ISR_multiplicity,h_jet_multiplicity,C,2,"Tracks","Num. Jets / Total");
	  C->Write();
	  C->Close();
		
	  salida_str = head_folder_results + "Delta_R_-_Jet_size";
	  salida = (Char_t *) salida_str.c_str();
	  C = new TCanvas(salida,"Delta R - Jet Size",1280,720);
	  Present(h_ISR_DeltaR,h_jet_DeltaR,C,1,"Delta_R","Num. Jets / Total");
	  C->Write();
	  C->Close();

	  // Correlated variables
	  salida_str = head_folder_results + "Cor_Delta_PT_Jet";
	  salida = (Char_t *) salida_str.c_str();
	  C = new TCanvas(salida,"Delta PT jet",1280,720);
	  Present(h_ISR_Delta_PT,h_jet_Delta_PT,C,2,"PT [GeV]","Num. Jets / Total");
	  C->Write();
	  C->Close();
		
	  salida_str = head_folder_results + "Cor_PT_proportion";
	  salida = (Char_t *) salida_str.c_str();
	  C = new TCanvas(salida,"PT proportion",1280,720);
	  Present(h_ISR_PT_HT,h_jet_PT_HT,C,2,"PT/HT","Num. Jets / Total");
	  C->Write();
	  C->Close();
		
	  salida_str = head_folder_results + "Cor_Delta_Eta_Average";
	  salida = (Char_t *) salida_str.c_str();
	  C = new TCanvas(salida,"Delta Eta Average",1280,720);
	  Present(h_ISR_Delta_Eta,h_jet_Delta_Eta,C,2,"Dh","Num. Jets / Total",122);
	  C->Write();
	  C->Close();
		
	  salida_str = head_folder_results + "Cor_Delta_Phi_Jet_MET_other_jets";
	  salida = (Char_t *) salida_str.c_str();
	  C = new TCanvas(salida,"Delta Phi - Jet MET - other jets",1280,720);
	  Present(h_ISR_DPhi_MET_other,h_jet_DPhi_MET_other,C,2,"Df","Num. Jets / Total",122);
	  C->Write();
	  C->Close();
		
	  salida_str = head_folder_results + "Cor_PT_over_<PT_other>";
	  salida = (Char_t *) salida_str.c_str();
	  C = new TCanvas(salida,"PT/<PT_other>",1280,720);
	  Present(h_ISR_PT_over_PT_others,h_jet_PT_over_PT_others,C,2,"PT/<PT>","Num. Jets / Total");
	  C->Write();
	  C->Close();

	  salida_str = head_folder_results + "Cor_Eta_over_<Eta_other>";
	  salida = (Char_t *) salida_str.c_str();
	  C = new TCanvas(salida,"Eta/<Eta_other>",1280,720);
	  Present(h_ISR_Eta_over_Eta_others,h_jet_Eta_over_Eta_others,C,3,"h/<h>","Num. Jets / Total",122);
	  C->Write();
	  C->Close();
	  
	  salida_str = head_folder_results + "Cor_Delta_Phi_over_<Delta_Phi_other>";
	  salida = (Char_t *) salida_str.c_str();
	  C = new TCanvas(salida,"Delta_Phi/<Delta_Phi_other>",1280,720);
	  Present(h_ISR_DPhi_over_Phi_others,h_jet_DPhi_over_Phi_others,C,3,"Df/<Df>","Num. Jets / Total",122);
	  C->Write();
	  C->Close();

	  // Comparison with the leading Jet
	  salida_str = head_folder_results + "Leading_Delta_PT";
	  salida = (Char_t *) salida_str.c_str();
	  C = new TCanvas(salida,"Delta PT: PT_leading-PT",1280,720);
	  Present(h_ISR_Delta_PT_leading,h_jet_Delta_PT_leading,C,2,"(PT_leading - PT)","Num. Jets / Total");
	  C->Write();
	  C->Close();

	  salida_str = head_folder_results + "Leading_Delta_Eta";
	  salida = (Char_t *) salida_str.c_str();
	  C = new TCanvas(salida,"Delta Eta: |Eta-Eta_leading|",1280,720);
	  Present(h_ISR_Delta_Eta_leading,h_jet_Delta_Eta_leading,C,2,"|Eta - Eta_leading|","Num. Jets / Total");
	  C->Write();
	  C->Close();
	  
	}
	
	hfile->Close();
	
	hfile_name_str = head_folder_results + "histos_Stops_200_N1_20-2.root";
	//hfile_name_str = head_folder_results + "histos_Stops_200_N1_35-2.root";
	//hfile_name_str = head_folder_results + "histos_Stops_300_N1_135-2.root";
	//hfile_name_str = head_folder_results + "histos_Stops_300_N1_120-2.root";
	//hfile_name_str = head_folder_results + "histos_Stops_400_N1_220-2.root";
	//hfile_name_str = head_folder_results + "histos_Stops_400_N1_235-2.root";
	//hfile_name_str = head_folder_results + "histos_Stops_500_N1_320-2.root";
	//hfile_name_str = head_folder_results + "histos-2.root";
	//hfile_name_str = head_folder_results + "histos_Top_semileptonico2.root";
	//hfile_name_str = head_folder_results + "histos_Top_Leptonico2.root";
	
	
	hfile_name = (Char_t *) hfile_name_str.c_str();

	TFile* hfile2 = new TFile(hfile_name, "RECREATE");
	h_ISR_PT_comp->Write();
	h_ISR_Eta_comp->Write();
	h_ISR_DPhi_MET_comp->Write();

	hist_ISR_PT->Write();
	hist_ISR_Abs_Eta->Write();
	hist_ISR_DPhi_MET->Write();
	hist_ISR_PT_ratio->Write();
	hist_ISR_Delta_Eta->Write();
	hist_ISR_DPhi_MET_other->Write();
	hist_ISR_Delta_PT_leading->Write();
	hist_ISR_Delta_Eta_leading->Write();

	hist_jet_PT->Write();
	hist_jet_Abs_Eta->Write();
	hist_jet_DPhi_MET->Write();
	hist_jet_PT_ratio->Write();
	hist_jet_Delta_Eta->Write();
	hist_jet_DPhi_MET_other->Write();
	hist_jet_Delta_PT_leading->Write();
	hist_jet_Delta_Eta_leading->Write();

	hfile2->Close();

	return 0;
}
