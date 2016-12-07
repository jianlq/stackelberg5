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

			CGraph *G = new CGraph("inputFile//testgraph.txt");
			vector<demand> reqbase;//background
			for(int i = 0; i < BGNUM; i++){
				int s = rand()%G->n, t;
				do{
					t = rand()%G->n;
				}while( s == t || G->canNotReach(s,t));
				reqbase.push_back(demand(s, t, rand()%(MAXDEMAND)+1));
			}

			/////dictor
			vector<demand> flow_OD;
			for(int i=0;i<overlay0->req.size();i++)
				flow_OD.push_back(overlay0->req[i]);
			for(int i=0;i<overlay1->req.size();i++)
				flow_OD.push_back(overlay0->req[i]);
			for(int i=0;i<overlay2->req.size();i++)
				flow_OD.push_back(overlay0->req[i]);
			for(int i=0;i<reqbase.size();i++)
				flow_OD.push_back(reqbase[i]);

			double eedic = 0;
			eedic = LBdictor(G,flow_OD);
			FILE *out = fopen("outputFile//lbdel.csv", "a");
			fprintf(out,"LB,%f,%f\n",eedic,G->delay);
			vector<double> ordic(overlay_num,0);
			
			ordic[0] = ORdictor(G,flow_OD,0,overlay0->req.size());
			fprintf(out,"OR0,%f,%f\n",G->mlu,ordic[0]);
			
			ordic[1] = ORdictor(G,flow_OD,overlay0->req.size(),overlay0->req.size()+overlay1->req.size());
			fprintf(out,"OR1,%f,%f\n",G->mlu,ordic[1]);
			
			ordic[2] = ORdictor(G,flow_OD,overlay0->req.size()+overlay1->req.size(),flow_OD.size()-reqbase.size());
			fprintf(out,"OR2,%f,%f\n",G->mlu,ordic[2]);
			fclose(out);


			//////初始分配
			overlay0->LP();
			overlay1->LP();	
			overlay2->LP();

			vector<vector<demand>> ALLreq;
			ALLreq.push_back(overlay0->getTrafficMatrix());
			ALLreq.push_back(overlay1->getTrafficMatrix());
			ALLreq.push_back(overlay2->getTrafficMatrix());

			vector<vector<demand>> overlayReq = ALLreq; ////before push_back background traffic ,backup overlay req for selfish NE

			ALLreq.push_back(reqbase);

			// initial flow_mark
			vector<int> overlayODnum;
			overlayODnum.push_back(overlay0->getEdgeNum());
			overlayODnum.push_back(overlay1->getEdgeNum());
			overlayODnum.push_back(overlay2->getEdgeNum());
			G->initMark(reqbase.size(),overlayODnum);
			
			//////relaxation for calculate NE 
			TE(G,reqbase,overlayReq);
			
			overlay0->getDelay(G->get_To_overlay(0));
			overlay1->getDelay(G->get_To_overlay(1));
			overlay2->getDelay(G->get_To_overlay(2));

			//////relaxation	
			vector<vector<demand>> updatereq;
			for(int i = 0;i < LOOP;i++){
				beta =(double)1/(i+1);
				updatereq.clear();

				overlay0->LP();
				overlay1->LP();	
				overlay2->LP();

				updatereq.push_back( overlay0->getTrafficMatrix());
				updatereq.push_back( overlay1->getTrafficMatrix());
				updatereq.push_back( overlay2->getTrafficMatrix());
		
				updatereq.push_back(reqbase);//background traffic

				G->CalculateDelay(&updatereq);
			
				overlay0->updateDelay(G->get_To_overlay(0));
				overlay1->updateDelay(G->get_To_overlay(1));
				overlay2->updateDelay(G->get_To_overlay(2));		
			}
			// util
			double util = G->CalculateMLU(&updatereq);


			G->initMark(reqbase.size(),overlayODnum);//mark置位
			//overlay延时置位  

			////Stackelberg SNE
			vector<overlay*> MutilOverlay;
			MutilOverlay.push_back(overlay0);
			MutilOverlay.push_back(overlay1);
			MutilOverlay.push_back(overlay2);	

			for(int con = 0;con < CONSIDER.size();con++)
			{
				evoluPopu popu(50, G->m, G, &ALLreq,MutilOverlay);
				evoluDiv ret = popu.evolution();
				cout << " WS "<<ret.ability<<endl;
				out = fopen("outputFile//lbdel.csv", "a");
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