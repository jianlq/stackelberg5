#include"LB.h"
#include"evolution.h"
#include"overlay.h"

int overlay_num = 3;
int LOOP = 100;
double alpha = 1;
double beta = 0.3;

int main(){
	srand((unsigned)time(NULL));
	int Time = 1;
	int CASEnum= 30;	

	vector<double>CONSIDER;
	int CON_VALUE = 10;
	double c =0;
	for(int i = 0; i <= CON_VALUE;i++){
		CONSIDER.push_back(c);
		c += 0.1;
	}

	for(int i =0;i<Time;i++){
		int conN = CONSIDER.size();
		vector<double> smlu(conN,0) ;
		vector<vector<double>> sdel(overlay_num);
		for(int i = 0;i<sdel.size();i++)
			sdel[i].resize(conN,0);
		vector<int> sucCaseSNE (conN, 0) ;

		double mlu = 0;
		vector<double> mludel(overlay_num, 0);

		int sucCaseLB = 0;

		vector<double> or(overlay_num,0);
		vector<double> ormlu(overlay_num,0);
		vector<int> sucCaseOR(overlay_num,0);

		double NEmlu = 0;
		vector<double> NEdel(overlay_num, 0);
		int sucCaseNE = 0;

		double NEmlu2 = 0;
		vector<double> NEdel2(overlay_num, 0);
		int sucCaseNE2 = 0;

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
			vector<demand> overlayreq = flow_OD;
			for(int i=0;i<reqbase.size();i++)
				flow_OD.push_back(reqbase[i]);

			FILE *out = fopen("outputFile//lbdel.csv", "a");
			fprintf(out,"\n\n\n casenum,%d \n",casenum);
			int end1 = overlay0->edge_flow.size();
			int end2 = overlay0->edge_flow.size()+overlay1->edge_flow.size();
			int end3 = flow_OD.size()-reqbase.size();
			vector<double> del = LBdictor(G,flow_OD,end1,end2,end3);
			if( G->mlu +1e-5 < INF){
				sucCaseLB++;
				mlu += G->mlu;
				mludel[0] += del[0];
				mludel[1] += del[1];
				mludel[2] += del[2];
				fprintf(out,"TE,,%f,%f,%f,%f\n",G->mlu,del[0],del[1],del[2]);
			}
			else
				break;

			/////// overlay0 optimal
		    vector<double> ordic(overlay_num,0);	
			vector<demand> background;
			for(int i=0;i<overlay1->edge_flow.size();i++)
				background.push_back(overlay1->edge_flow[i]);
			for(int i=0;i<overlay2->edge_flow.size();i++)
				background.push_back(overlay2->edge_flow[i]);
			for(int i=0;i<reqbase.size();i++)
				background.push_back(reqbase[i]);
			ordic[0] = ORdictor(G,background,overlay0->edge_flow);
			if( ordic[0] + 1e-5 <INF){
				sucCaseOR[0]++;
				or[0] += ordic[0];
				ormlu[0] += G->mlu;
				fprintf(out,"OR0,,%f,%f\n",G->mlu,ordic[0]);
			}
			else
				break;

			/////// overlay1 optimal
			background.clear();
			for(int i=0;i<overlay0->edge_flow.size();i++)
				background.push_back(overlay0->edge_flow[i]);
			for(int i=0;i<overlay2->edge_flow.size();i++)
				background.push_back(overlay2->edge_flow[i]);
			for(int i=0;i<reqbase.size();i++)
				background.push_back(reqbase[i]);
			ordic[1] = ORdictor(G,background,overlay1->edge_flow);
			if( ordic[1] + 1e-5 < INF ){
				sucCaseOR[1]++;
				or[1] += ordic[1];
				ormlu[1] += G->mlu;
				fprintf(out,"OR1,,%f,,%f\n",G->mlu,ordic[1]);
			}
			else
				break;
		
			
			/////// overlay2 optimal
			background.clear();
			for(int i=0;i<overlay0->edge_flow.size();i++)
				background.push_back(overlay0->edge_flow[i]);
			for(int i=0;i<overlay1->edge_flow.size();i++)
				background.push_back(overlay1->edge_flow[i]);
			for(int i=0;i<reqbase.size();i++)
				background.push_back(reqbase[i]);
			ordic[2] = ORdictor(G,background,overlay2->edge_flow);
			if( ordic[2] + 1e-5 <INF){
				sucCaseOR[2]++;
				or[2] += ordic[2];
				ormlu[2] += G->mlu;
				fprintf(out,"OR2,,%f,,,%f\n",G->mlu,ordic[2]);
			}
			else
				break;


			double or = ORdictor(G,background,overlayreq);
			fprintf(out,"OR,,%f,,,%f\n",G->mlu,or);

			// initial flow_mark
			vector<int> overlayODnum;
			overlayODnum.push_back(overlay0->getEdgeNum());
			overlayODnum.push_back(overlay1->getEdgeNum());
			overlayODnum.push_back(overlay2->getEdgeNum());
			G->initMark(reqbase.size(),overlayODnum);
			
			bool can = TE(G,reqbase,overlayReq);////////relaxation for calculate NE 
			G->CalculateDelay(&ALLreq);
			if( can ){
				overlay0->getDelay(G->get_To_overlay(0));
				overlay1->getDelay(G->get_To_overlay(1));
				overlay2->getDelay(G->get_To_overlay(2));

				//////relaxation	
				vector<vector<demand>> updatereq;
				for(int i = 0;i < 100;i++){
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
				fprintf(out,"NE,,%f,%f,%f,%f\n",util,overlay0->getCost(),overlay1->getCost(),overlay2->getCost());
				
				
				NEmlu += util;
				NEdel[0] += overlay0->getCost();
				NEdel[1] += overlay1->getCost();
				NEdel[2] += overlay2->getCost();
				sucCaseNE++;
			}
			else
				break;

			G->clearMark();//mark置位 
			////初始TE用的是规划 以网络目标最小，得到流分配
			double ne = TEcplex(G,reqbase,overlayReq); ////////relaxation for calculate NE 		
			if( ne + 1e-5 <INF ){
				overlay0->getDelay(G->get_To_overlay(0));
				overlay1->getDelay(G->get_To_overlay(1));
				overlay2->getDelay(G->get_To_overlay(2));

				//////relaxation	
				vector<vector<demand>> updatereq;
				for(int i = 0;i < 100;i++){
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
				fprintf(out,"NE2,,%f,%f,%f,%f\n",util,overlay0->getCost(),overlay1->getCost(),overlay2->getCost());
				fclose(out);
				
				NEmlu2 += util;
				NEdel2[0] += overlay0->getCost();
				NEdel2[1] += overlay1->getCost();
				NEdel2[2] += overlay2->getCost();
				sucCaseNE2++;
			}
			else
				break;


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
			double sumcost = NEdel[0]+NEdel[1]+NEdel[2];
			for(int con = 0;con < CONSIDER.size();con++)
			{
				evoluPopu popu(50, G->m, G, &ALLreq,MutilOverlay,CONSIDER[con]*sumcost);
				evoluDiv ret = popu.evolution();
				if(ret.ability + 1e-5 < INF ){
					smlu[con] += ret.mlu;
					sdel[0][con] += overlay0->getCost();
					sdel[1][con] += overlay1->getCost();
					sdel[2][con] += overlay2->getCost();
					sucCaseSNE[con]++;
					out = fopen("outputFile//lbdel.csv", "a");
					fprintf(out,"SNE,%f,%f,%f,%f,%f\n",CONSIDER[con],ret.ability,overlay0->getCost(),overlay1->getCost(),overlay2->getCost());
					fclose(out);
				}
			}
			delete G;
			delete overlay0;
			delete overlay1;
			delete overlay2;
		} // end of CaseNum for

		FILE *res = fopen("outputFile//avgresult.csv", "a");		
		fprintf(res,"\n\n%d case average\n",CASEnum);
		fprintf(res,",CONSIDER,sucCase,mlu,overlay0 delay,overlay1 delay,overlay2 delay\n");

		fprintf(res, "LB,,%d,%lf,%lf,%lf,%lf\n",sucCaseLB,
			mlu/sucCaseLB,
			mludel[0]/sucCaseLB,
			mludel[1]/sucCaseLB,
			mludel[2]/sucCaseLB);

		fprintf(res, "OR0,,%d,%lf,%lf\n",sucCaseOR[0],
			ormlu[0]/sucCaseOR[0],
			or[0]/sucCaseOR[0]);

		fprintf(res, "OR1,,%d,%lf,,%lf\n",sucCaseOR[1],
			ormlu[1]/sucCaseOR[1],
			or[1]/sucCaseOR[1]);

		fprintf(res, "OR2,,%d,%lf,,,%lf\n",sucCaseOR[2],
			ormlu[2]/sucCaseOR[2],
			or[2]/sucCaseOR[2]);

		fprintf(res, "NE,,%d,%lf,%lf,%lf,%lf\n",sucCaseNE,
			NEmlu/sucCaseNE,
			NEdel[0]/sucCaseNE,
			NEdel[1]/sucCaseNE,
			NEdel[2]/sucCaseNE );

		fprintf(res, "NE2,,%d,%lf,%lf,%lf,%lf\n",sucCaseNE2,
			NEmlu2/sucCaseNE2,
			NEdel2[0]/sucCaseNE2,
			NEdel[1]/sucCaseNE2,
			NEdel2[2]/sucCaseNE2 );
		
		for(unsigned int con = 0;con < CONSIDER.size();con++){
			fprintf(res, "SNE,%lf,%d,%lf,%lf,%lf,%lf\n",CONSIDER[con],sucCaseSNE[con],
				smlu[con]/sucCaseSNE[con],
				sdel[0][con]/sucCaseSNE[con],
			    sdel[1][con]/sucCaseSNE[con],
			    sdel[2][con]/sucCaseSNE[con]); 
		}
		fclose(res);
	} //// end of Time for
	system("pause");
	return 0;	
}