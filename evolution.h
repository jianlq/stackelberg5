﻿#ifndef EVOLUTION
#define EVOLUTION
#include "Common.h"
#include"CGraph.h"
#include"overlay.h"

extern int LOOP;
class evoluDiv{
	private:
		static const int MUT = 10;
		static const int CUL = 30;
		static const int HOR = 10;
		CGraph *G;
		vector<vector<demand>> *dem;
		vector<overlay*> Overlay;
		double consider;
	public:
		vector<double> x;
		double mlu;
		double cost;
		double ability;
		evoluDiv() {;}

		void init(){
			ability = INF;
			mlu = INF;
			cost = INF;
		}
		evoluDiv(int m, CGraph *g, vector<vector<demand>> *d,vector<overlay*>&ov,double con){
			init();
			x.resize(m);
			G = g;
			dem = d;
			Overlay = ov;
			consider = con;
			randNature();	
		}

		evoluDiv(vector<double> &tx, CGraph *g, vector<vector<demand>> *d,vector<overlay*> &ov,double con){
			init();
			x.clear();
			G = g;
			dem = d;
			Overlay = ov;
			consider = con;
			for(int i = 0; i < tx.size(); i++)
				x.push_back(tx[i]);
		}

		evoluDiv mate(evoluDiv other){ 
			vector<double> nx;
			for(int i = 0; i < x.size(); i++){
				double r = 1.0 * rand() / (RAND_MAX);
				nx.push_back(r * x[i] + (1 - r) * other.x[i]);
			}
			return evoluDiv(nx, G, dem,Overlay,consider);
		}

		bool TE(){
			int can = 0;
			int sumreq = 0;
			int num = (*dem).size(); //OR1,OR2,OR3,background
			for(int i = 0; i < num ; i++)
				sumreq += (*dem)[i].size();	

			for(int d = 0; d < (*dem)[num-1].size();d++){
				double dis = G->dijkstraWeight(d,(*dem)[num-1][d].org, (*dem)[num-1][d].des, (*dem)[num-1][d].flow);
				if( dis+1e-5 <=INF){
					for( int k = 0;k< G->reqPathID[d].size();k++){
						G->background_mark[d][G->reqPathID[d][k]] = 1;
					}
					can++;
				}
			}

			for(int i = 0; i < num-1; i++){  //// i:Overlay No
				for(int d = 0; d < (*dem)[i].size();d++){
					double dis = G->dijkstraWeight(d,(*dem)[i][d].org, (*dem)[i][d].des, (*dem)[i][d].flow);
					if( dis+1e-5 < INF){
						for( int k = 0; k < G->reqPathID[d].size();k++){
							G->overlay_mark[i][d][G->reqPathID[d][k]] = 1;
						}
						can++;
					}
				}
			}

			if( can < sumreq ){
				return false;
			}
			return true;
		}

		////ability
		void calAbility(){  
			mlu = 0;
			ability = 0;
			G->clearOcc();
			for(int i = 0; i < G->m; i++)
				G->Link[i]->dist = x[i];

			// relaxation algorithm
			if(!TE()){
				ability = INF;
				mlu = INF;
				return ;
			}
			//underlay calculate delay for overlay
			G->CalculateDelay(dem);
			for(int i = 0; i < Overlay.size(); i++){
				(*Overlay[i]).getDelay(G->get_To_overlay(i));
			}

			//////	relaxation	for calculate NE	
			vector<vector<demand>> updatereq;
			for(int i = 0;i < LOOP;i++){
				beta =(double)1/(i+1);
				updatereq.clear();
				for(int i = 0; i < Overlay.size(); i++){
					Overlay[i]->LP();
					updatereq.push_back( Overlay[i]->getTrafficMatrix());
				}
				updatereq.push_back((*dem)[(*dem).size()-1]);//background traffic
				G->CalculateDelay(&updatereq);
				for(int i = 0; i < Overlay.size(); i++){
					Overlay[i]->updateDelay(G->get_To_overlay(i));
				}
			}

			// util
			double util = G->CalculateMLU(&updatereq);
			mlu = util;
			cout << " mlu: \t" << util <<"overlay0:"<< Overlay[0]->getCost() << " \toverlay1:" << Overlay[1]->getCost() << " \toverlay2:" << Overlay[2]->getCost()<<endl;
			
			double sumdelay = Overlay[0]->getCost() + Overlay[1]->getCost() +Overlay[2]->getCost() ;
			ability = util + sumdelay/this->consider;
			//ability = util;
			//ability = 1.0/sumdelay;
			cost = sumdelay;

			G->clearMark();
		}

		void randNature(){
			for(int i = 0; i < x.size(); i++){
				x[i] = 1.0 * HOR * rand() / RAND_MAX;
				if(rand()%2)
					x[i] = -x[i];
				x[i] += G->Link[i]->dist;
				x[i] = max(0.0, x[i]);
				x[i] = min(1.0*MAXWEIGHT, x[i]);
			}
		}

		void mutation(){  //变异
			for(int i = 0; i < x.size(); i++){
				x[i] += -(MAXWEIGHT * MUT/100.0) + 2 * (MAXWEIGHT * MUT/100.0) * rand() / RAND_MAX;
				x[i] = max(0.0, x[i]);
				x[i] = min(1.0 * MAXWEIGHT, x[i]);
			}
		}
		void culture(evoluDiv hero){
			for(int i = 0; i < x.size(); i++){
				double r = (CUL/100.0) * rand() / RAND_MAX;
				x[i] = (1 - r) * x[i] + r * hero.x[i];
			}
		}
};

bool evoluCmp(evoluDiv a, evoluDiv b){
	return a.ability < b.ability;
}

class evoluPopu{
	private:
		static const int YEAR = 60;
		static const int YEARCUL = 30;
		static const int NOHERO = 20;
		vector<evoluDiv> popu;
		CGraph *G;
		vector<vector<demand>> *dem;
		FILE *herofile;
	public:
		evoluDiv hero;
		// n 个体数，m：每个个体对应的解
		evoluPopu(int n, int m, CGraph *g, vector<vector<demand>> *d,vector<overlay*> &mutil,double con){
			popu.clear();
			G = g;
			dem = d;
			for(int i = 0; i < n; i++){
				evoluDiv divi(m, G, dem,mutil,con);
				popu.push_back(divi);
			}
			hero = evoluDiv(m, G, dem,mutil,con);
			herofile = fopen("outputFile//hero.csv","a");
		}
		evoluDiv evolution(){
			double costlog = INF;
			int index = 0;

			int nohero = 0;
			fprintf(herofile,"Start:\n ");
			for(int i = 0; i < hero.x.size(); i++)
				hero.x[i] = G->Link[i]->dist;
			hero.calAbility();
			fprintf(herofile, "%f,%f\n", hero.mlu,hero.cost);

			for(int i = 0; i < popu.size(); i++){
				popu[i].calAbility();
				if(popu[i].cost < costlog){
					costlog = popu[i].cost;
					index = i;
				}
			}
			fprintf(herofile, "%f,%f\n", popu[index].mlu,costlog);

			sort(popu.begin(), popu.end(), evoluCmp);
			
			for(int curYear = 1; curYear <= YEAR; curYear++){
				int n = popu.size(), getMore = 0;
				vector<evoluDiv> sons;
				for(int i = 0; i+1 < n; i+=2){
					sons.push_back(popu[i].mate(popu[i+1]));
					sons.push_back(popu[i+1].mate(popu[i]));
					sons.push_back(popu[i].mate(popu[i+1]));
				}
				int m = sons.size();
				for(int i = 0; i < m; i++){
					sons[i].mutation();
					if(curYear > YEARCUL)
						sons[i].culture(hero);
					sons[i].calAbility();
				}
				sort(sons.begin(), sons.end(), evoluCmp);
				popu.clear();
				for(int i = 0; i < n; i++){
					popu.push_back(sons[i]);
					if(sons[i].ability < hero.ability){
						hero = sons[i];
						getMore = 1;
					}

					if(sons[i].cost < costlog){
						costlog = sons[i].cost;
						index = i;
					}
					
				}
				fprintf(herofile, "Year %d:\n", curYear);
				fprintf(herofile,"%f,%f\n",sons[index].mlu,costlog);

				if(getMore){
					;
				}
				else 
					nohero++;
				if(nohero> NOHERO){
					break;
				}
			}
			fprintf(herofile,"end\n\n\n\n");
			fclose(herofile);
			return hero;
		}
};

#endif
