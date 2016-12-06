#include"LB.h"
#include"evolution.h"
#include"overlay.h"

int overlay_num = 3;
int LOOP = 50;
double alpha = 1;
double beta = 0.3;



int main(){
	srand((unsigned)time(NULL));
	int Time = 3;
	int CASEnum= 30;	

	vector<double>CONSIDER;
	int CON_VALUE = 25;
	double c =0;
	for(int i=0;i <= CON_VALUE;i++){
		CONSIDER.push_back(c);
		c += 0.2;
	}

	for(int i =0;i<Time;i++){

		FILE *out = fopen("outputFile//lbdel.csv", "a");

		int conN = CONSIDER.size();
		vector<double> cmpmlu(conN,0) ;
		vector<double> cmpdelay(conN,0);
		vector<int> successCase (conN, 0) ;

		for(int casenum = 0; casenum < CASEnum; casenum++){

			overlay *overlay0 = new overlay(0,"inputFile//overlay0.txt");
			overlay *overlay1 = new overlay(1,"inputFile//overlay1.txt");
			overlay *overlay2 = new overlay(2,"inputFile//overlay2.txt");

			overlay0->get_OD_infor("inputFile//overlay0_OD.txt");
			overlay1->get_OD_infor("inputFile//overlay1_OD.txt");
			overlay2->get_OD_infor("inputFile//overlay2_OD.txt");

			overlay0->LP();
			overlay1->LP();
			overlay2->LP();

			vector<vector<demand>> ALLreq;
			ALLreq.push_back(overlay0->getTrafficMatrix());
			ALLreq.push_back(overlay1->getTrafficMatrix());
			ALLreq.push_back(overlay2->getTrafficMatrix());

			CGraph *G = new CGraph("inputFile//graphATT.txt");


			vector<demand> reqbase;//background
			for(int i = 0; i < BGNUM; i++){
				int s = rand()%G->n, t;
				do{
					t = rand()%G->n;
				}while( s == t || G->canNotReach(s,t));
				reqbase.push_back(demand(s, t, rand()%(MAXDEMAND)+1));
			}
			ALLreq.push_back(reqbase);

			// initial flowmark
			vector<int> overlayODnum;
			overlayODnum.push_back(overlay0->getEdgeNum());
			overlayODnum.push_back(overlay1->getEdgeNum());
			overlayODnum.push_back(overlay2->getEdgeNum());
			G->initMark(reqbase.size(),overlayODnum);


			vector<demand> req;
			for(int i=0;i<ALLreq.size();i++)
				for(int d =0;d <ALLreq[i].size();d++)
					req.push_back(ALLreq[i][d]);

			double eedic = 0;
			eedic = LBdictor(G,req);
			fprintf(out,"EE,%f,%f\n",eedic,G->delay);
			vector<double> ordic(overlay_num,0);
			
			ordic[0] = ORdictor(G,req,0,overlay0->req.size());
			fprintf(out,"OR0,%f,%f\n",G->mlu,ordic[0]);
			
			ordic[1] = ORdictor(G,req,overlay0->req.size(),overlay0->req.size()+overlay1->req.size());
			fprintf(out,"OR1,%f,%f\n",G->mlu,ordic[1]);
			
			ordic[2] = ORdictor(G,req,overlay0->req.size()+overlay1->req.size(),req.size()-reqbase.size());
			fprintf(out,"OR2,%f,%f\n",G->mlu,ordic[2]);
		
			vector<overlay*> MutilOverlay;
			MutilOverlay.push_back(overlay0);
			MutilOverlay.push_back(overlay1);
			MutilOverlay.push_back(overlay2);		
			
		
			for(int con = 0;con < CONSIDER.size();con++){

				evoluPopu popu(50, G->m, G, &ALLreq,MutilOverlay);
				evoluDiv ret = popu.evolution();
				cout << " WS "<<ret.ability<<endl;
				fprintf(out,"WS,%f,%f,%f,%f\n",ret.ability,overlay0->getCost(),overlay1->getCost(),overlay2->getCost());
		    	fclose(out);
				exit(1);

			}

			delete G;
			delete overlay0;
			delete overlay1;
			delete overlay2;

		}
	/*	FILE *res = fopen("outputFile//result.csv", "a");

		fclose(res);*/

	}
	system("pause");
	return 0;	
}