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

			CGraph *G = new CGraph("inputFile//ATT.txt");
			vector<demand> reqbase;//background
			for(int i = 0; i < BGNUM; i++){
				int s = rand()%G->n, t;
				do{
					t = rand()%G->n;
				}while( s == t || G->canNotReach(s,t));
				reqbase.push_back(demand(s, t, rand()%(MAXDEMAND)+1));
			}

			FILE *out = fopen("outputFile//lbdel.csv", "a");
			///////each overlay selfish规划
		     vector<double> ordic(overlay_num,0);	
			vector<demand> background0;
			for(int i=0;i<overlay1->req.size();i++)
				background0.push_back(overlay1->req[i]);
			for(int i=0;i<overlay2->req.size();i++)
				background0.push_back(overlay2->req[i]);
			for(int i=0;i<reqbase.size();i++)
				background0.push_back(reqbase[i]);
			ordic[0] = ORdictor(G,background0,overlay0->req);
			fprintf(out,"OR0,%f,%f\n",G->mlu,ordic[0]);
			

			vector<demand> background1;
			for(int i=0;i<overlay0->req.size();i++)
				background1.push_back(overlay0->req[i]);
			for(int i=0;i<overlay2->req.size();i++)
				background1.push_back(overlay2->req[i]);
			for(int i=0;i<reqbase.size();i++)
				background1.push_back(reqbase[i]);
			ordic[1] = ORdictor(G,background1,overlay1->req);
			fprintf(out,"OR1,%f,%f\n",G->mlu,ordic[1]);
			
			vector<demand> background2;
			for(int i=0;i<overlay0->req.size();i++)
				background2.push_back(overlay0->req[i]);
			for(int i=0;i<overlay1->req.size();i++)
				background2.push_back(overlay1->req[i]);
			for(int i=0;i<reqbase.size();i++)
				background2.push_back(reqbase[i]);
			ordic[2] = ORdictor(G,background2,overlay2->req);
			fprintf(out,"OR2,%f,%f\n",G->mlu,ordic[2]);
			//exit(1);


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

			/////LBdictor flow_OD
			vector<demand> flow_OD;
			for(int i=0;i<overlay0->edge_flow.size();i++)
				flow_OD.push_back(overlay0->edge_flow[i]);
			for(int i=0;i<overlay1->edge_flow.size();i++)
				flow_OD.push_back(overlay1->edge_flow[i]);
			for(int i=0;i<overlay2->edge_flow.size();i++)
				flow_OD.push_back(overlay2->edge_flow[i]);
			for(int i=0;i<reqbase.size();i++)
				flow_OD.push_back(reqbase[i]);

			double eedic = 0;
			int end1 = overlay0->edge_flow.size(),end2=overlay0->edge_flow.size()+overlay1->edge_flow.size();
			int end3 = flow_OD.size()-reqbase.size();
			eedic = LBdictor(G,flow_OD);
			fprintf(out,"\n\nLB,%f,%f\n",eedic,G->delay);

		
			ordic[0] = ORdictor(G,flow_OD,0,end1);
			fprintf(out,"OR0,%f,%f\n",G->mlu,ordic[0]);
			
			ordic[1] = ORdictor(G,flow_OD,end1,end2);
			fprintf(out,"OR1,%f,%f\n",G->mlu,ordic[1]);
			
			ordic[2] = ORdictor(G,flow_OD,end2,end3);
			fprintf(out,"OR2,%f,%f\n",G->mlu,ordic[2]);

			// initial flow_mark
			vector<int> overlayODnum;
			overlayODnum.push_back(overlay0->getEdgeNum());
			overlayODnum.push_back(overlay1->getEdgeNum());
			overlayODnum.push_back(overlay2->getEdgeNum());
			G->initMark(reqbase.size(),overlayODnum);

			//////relaxation for calculate NE 	
			 TE(G,reqbase,overlayReq);
			 fprintf(out,"TE,%f\n",G->mlu);

			overlay0->getDelay(G->get_To_overlay(0));
			overlay1->getDelay(G->get_To_overlay(1));
			overlay2->getDelay(G->get_To_overlay(2));

			//////relaxation	
			vector<vector<demand>> updatereq;
			for(int i = 0;i < 200;i++){
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
			cout << updatereq.size();
			// util
			double util = G->CalculateMLU(&updatereq);
			cout << "NE "<<" util "<<util<<endl;
			fprintf(out,"NE,%f,%f,%f,%f\n",util,overlay0->getCost(),overlay1->getCost(),overlay2->getCost());
			fclose(out);
			exit(1);

			G->clearMark();//mark置位
			
			//overlay延时置位
			overlay0->clearDelay();
			overlay1->clearDelay();
			overlay2->clearDelay();

			////Stackelberg SNE
			vector<overlay*> MutilOverlay;
			MutilOverlay.push_back(overlay0);
			MutilOverlay.push_back(overlay1);
			MutilOverlay.push_back(overlay2);	

			for(int con = 0;con < CONSIDER.size();con++)
			{
				evoluPopu popu(100, G->m, G, &ALLreq,MutilOverlay);
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