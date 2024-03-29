# ISR_uniandes

## **Instructions to copy the repository**
1. Create an account in github
2. Set up your user details in your laptop/dzero:
	* git config --global user.name "Your github Name here"
	* git config --global user.email "your_email@domainname.com"
3. Set up you ssh key:
	* Follow the instructions: https://help.github.com/articles/generating-ssh-keys/#platform-linux
4. Create a folder with the name of the repository:
	* mkdir ISR_uniandes
5. Make a local copy of the repository:
	* git clone git@github.com:jotadram6/ISR_uniandes.git ISR_uniandes

## **To edit this README file**
Follow the instructions here: https://guides.github.com/features/mastering-markdown/

## **Minimal instructions**
A git repository is based on retrivieng (pull) and sending (push) from and for the repository.

_Don't forget to **always pull before start working** and to **commit and push when you are done**!!!_

This repository is intended to share code, please don't submit any samples (root, lhe, hep files)!
* Get all the changes from other contributors: git pull origin master
* Add the new files you want to submit to the repository: git add _filename_
* Commit your changes: git commit -a -m "Your commit message"
	* Don't be shy with your message! Clear (long) messages allow a better collaboration. 
* And push! : git push origin master

In order to create and submit to your own branch, please follow these instructions:
* Create your branch: git branch branch_name
* Get the list of branches: git branch
* Switch to your branch: git checkout branch_name
* Pull and push from your branch:
	* git pull origin branch_name
	* git push origin branch_name

## **Running the code**
1) Se debe realizar el correspondiente ajuste de los directorios en donde van a quedar los .root en el archivo ISR_matching.cpp y en el config_file.txt para el 
ISR_matching, ISR_tagging y el ISR_jet_analysis.
2) Primero se compila el ISR_Matching con el comando "make compile_ROOT_Delphes" que permite compilar ROOT con DELPHES.
3) Ahora se corre el ejecutable con "./ISR_matching config_file.txt".
4) Se realiza el mismo procedimiento de compilar y correr el ejecutable para el ISR tagging ("make compile_ROOT_Delphes" y luego "./ISR_tagging config_file.txt")
 y el ISR jet Analysis ("make compile_ROOT_Delphes" y luego "./ISR_jet_analysis config_file.txt").
5) Es importante realizar las compilaciones de esta forma lineal, teniendo en cuenta que el ISR jet Analysis corre sobre los resultados arrojados por el algoritmo de 
matching.

## **Bibliography**

1. Final report from Andres Felipe Garcia:
	* https://github.com/andresfgarcia150/ISR_tagging_project/blob/master/Documentation/Final_report.pdf
2. David Krohn and Lisa Randal:
	* http://arxiv.org/pdf/1101.0810v1.pdf
3. Teruki:
	* PRD 92, 095009
	* PRD 90, 095002
4. http://arxiv.org/pdf/1506.07885v1.pdf
5. Lian-Tao:
	* http://arxiv.org/pdf/1506.00653v1.pdf

CMS-related:

1. XS's:
	* https://twiki.cern.ch/twiki/bin/view/LHCPhysics/SUSYCrossSections
	* https://twiki.cern.ch/twiki/bin/view/CMSPublic/WinoCn
	* http://pauli.uni-muenster.de/~akule_01/nllwiki/index.php/NLL-fast
2. https://twiki.cern.ch/twiki/bin/view/CMS/SMSTChiSlepSnuMadgraph8TeV
3. https://twiki.cern.ch/twiki/bin/view/CMS/SMST2DegenerateStopMadgraph8TeV
4. https://twiki.cern.ch/twiki/bin/viewauth/CMS/SMST2ccMadgraph8TeV
5. Analysis: https://cds.cern.ch/record/2010110/files/SUS-14-021-pas.pdf

