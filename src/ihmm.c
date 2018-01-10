
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "tldevel.h"

#include "ihmm.h"
#include <math.h>
#include <float.h>

struct pgas{
        int** particle_ancestry;
        float* x;
        float* y;
        float* particle_path_prob;
        float* tmp_particle_path_prob;
        float* previous_path_score;
        int* ancestors;
        float* ancestor_weight;
        int num_particles;
        int buff_size;
};

#define EMISSION_H 0.3

struct dp{
        int* state_used;
        int* ID;
        float** dp_matrix;
        int alloced_states;
        int alloced_len;
};

struct ihmm_sequences{
        char** seq;
        float** u;
        float* score;
        int** labels;
        int* len;
        int num_seq;
        int max_len;
};

static int select_random(float* vector,int num_elem,int* selected,float r);

static int iHmmSampleBeam(struct ihmm_sequences* iseq,struct iHMM_model* model);
static int iHmmHyperSample(struct iHMM_model* model, int iterations);

static int set_u(struct iHMM_model* model, struct ihmm_sequences* iseq, float* min_u);

static int dyn_prog(struct dp* dyn,struct iHMM_model* model, struct ihmm_sequences* iseq,int num);
static int remove_unused_states(struct dp* dyn,struct iHMM_model* model, struct ihmm_sequences* iseq);
static int print_hyper_parameters(struct iHMM_model* model);

static struct ihmm_sequences* init_ihmm_seq(char** sequences,int numseq);
static void free_ihmm_seq(struct ihmm_sequences* iseq);
static int set_labels_based_on_alignment(struct ihmm_sequences* iseq,char** aln,int aln_len,int num_seq_in_aln);

static int add_state_to_model(struct iHMM_model* model);

static struct dp* init_dp(int num_states, int len);
static int clear_dp(struct dp* dp);
static int resize_dp(struct dp* dp, int num_states,int len);
static int free_dp(struct dp* dp);

static int add_random_counts(struct iHMM_model* model, int count);
static int fill_counts(struct ihmm_sequences* iseq,struct iHMM_model* model);
static int sample_counts(struct iHMM_model* model);
static int clear_counts(struct iHMM_model* model);
static int trim_counts(struct iHMM_model* model);

static int print_prior(struct iHMM_model* model);
static int print_transistion_matrix(struct iHMM_model* model);
static int print_trans_count_matrix(struct iHMM_model* model);
static int print_emisson_matrix(struct iHMM_model* model);

static int pgas_sample(struct iHMM_model* model, struct pgas* pgas,struct ihmm_sequences* iseq,int num);

static int remove_unused_states_smc(struct iHMM_model* model, struct ihmm_sequences* iseq);

static struct pgas* init_pgas(int max_len, int malloced_states,int num_particles);
static int resize_pgas(struct pgas* pgas, int malloced_states);
static void free_pgas(struct pgas* pgas);

int run_make_ihmm(struct iHMM_model* model,char** sequences,int numseq)
{
        struct ihmm_sequences* iseq = NULL;
        //struct model* model = NULL;
	
        int aln_len = 0;
        int i,j;
        ASSERT(numseq !=0,"no sequences provided");
	
        RUNP(iseq =init_ihmm_seq(sequences,numseq));
	
        //sequences = kalign_align(sequences, numseq);
	
        aln_len = (int) strlen(sequences[0]);
	
        RUN(set_labels_based_on_alignment(iseq, sequences, aln_len,numseq));
	
        //testing
        /*rk_randomseed( &rndstate);
        //aln_len = 4;
        for(i = 0; i < numseq;i++){
        c =0;
        for(j = 0; j < aln_len;j++){
        switch (sequences[i][j]){
        case 'A':
        //	iseq->labels[i][c] = 0;
        //	fprintf(stdout,"%3d",iseq->labels[i][c]);
        //	c++;
        //	break;
        case 'C':
        //	iseq->labels[i][c] = 1;
        //	fprintf(stdout,"%3d",iseq->labels[i][c]);
        //	c++;
        //	break;
        case 'G':
        ///	iseq->labels[i][c] = 2;
        //	fprintf(stdout,"%3d",iseq->labels[i][c]);
        //	c++;
        //	break;
        case 'T':
	 
        iseq->labels[i][c] = j;
        //iseq->labels[i][c] =  (int) rk_interval(3-1, &rndstate);
        fprintf(stdout,"%3d",iseq->labels[i][c]);
        c++;
        break;
        default:
        fprintf(stdout,"   ");
        break;
        }
	 
        }
        fprintf(stdout,"\n");
        }
        */
        // */
        //RUNP(model = init_model(),"init hypers failed.");
	
        // set expected states to number of columns in the kalign alignment
        model->expected_K = aln_len;
	
        //for(i = 0; i < numseq;i++){
        //	for(j = 0; j < iseq->len[i];j++){
        //iseq->labels[i][j] = (int)rk_interval(model->expected_K -1,   &model->rndstate);
        //	}
        //}
        //
	
	
	
        //####SUPER critical - otherwise ihmm code looks for states that do not exist....
	
        RUN(iHmmSampleBeam(iseq, model));
	
	
	
	
        //print_transistion_matrix(model);
//	print_emisson_matrix(model);
	
	
	
        for(i = 0; i < 4;i++){
                model->back[i] = 0.0f;
        }
        for(i = 0; i < numseq;i++){
                for(j = 0; j < iseq->len[i];j++){
			
                        model->back[(int) iseq->seq[i][j]]++;
                }
        }
	
	
        float sum = 0.0;
	
        for(i = 0; i < 4;i++){
                sum += model->back[i];// 0.0f;
        }
        for(i = 0; i < 4;i++){
                model->back[i] /= sum;
                model->back[i] = prob2scaledprob(model->back[i]);
        }

	
	
        /*model->K = 1;
          model->transition[0][0] = 0.9;
          model->transition[0][1] = 0.1;
          model->transition[1][0] = 0.1;
          model->transition[1][1] = 0.9;
	 
          model->emission[0][0] = 0.8;
          model->emission[0][1] = 0.1;
          model->emission[0][2] = 0.1;
	 
	 
          model->emission[1][0] = 0.1;
          model->emission[1][1] = 0.1;
          model->emission[1][2] = 0.8;
	 
	 
	 
          iseq->seq[0][0] = 0;
          iseq->seq[0][1] = 0;
          iseq->seq[0][2] = 0;
          iseq->seq[0][3] = 1;
          iseq->seq[0][4] = 1;
          iseq->seq[0][5] = 1;
          iseq->seq[0][6] = 2;
          iseq->seq[0][7] = 2;
          iseq->seq[0][8] = 2;*/
	
	
        //smc_sampler(model, iseq,0);// iseq->len[0]); //??????
	
        /*c = 0;
          for(i = 0; i < numseq;i++){
          for(j = 0; j < aln_len;j++){
          fprintf(stdout,"  %c",sequences[i][j]);
          }
          fprintf(stdout,"\n");
          c = 0;
          for(j = 0; j < aln_len;j++){
          switch (sequences[i][j]){
          case 'A':
          case 'C':
          case 'G':
          case 'T':
					
					//iseq->labels[i][c] = j;
					fprintf(stdout,"%3d",iseq->labels[i][c]);
					c++;
					break;
          default:
					fprintf(stdout,"   ");
					break;
          }
          }
          fprintf(stdout,"\n");
		
          }
          fprintf(stdout,"\n");*/
	
        //free_model(model);
	
	
        free_ihmm_seq(iseq);
        return OK;
ERROR:
        free_ihmm_seq(iseq);
        return FAIL;
}


int particle_gibbs_with_ancestors_controller(struct iHMM_model* model,char** sequences,int numseq)
{
        struct ihmm_sequences* iseq = NULL;
        struct pgas* pgas = NULL;
        float sum = 0.0;
        
        int i,j;
	
        int iter;
	
        ASSERT(numseq !=0,"no sequences provided");
        init_logsum();
		
        RUNP(iseq =init_ihmm_seq(sequences,numseq));
	

        if(model->expected_K == 0){
                WARNING_MSG("Expected number of states is set to zero.");
                model->expected_K = iseq->max_len*5;
                LOG_MSG("Setting expected number of states to %d.",model->expected_K);
                            
        }

        
        /* Setting up model structure */
        
        RUN(start_iHMM_model(model, model->expected_K));

        /* Setting up structure do deal with particles  */
        RUNP(pgas = init_pgas(iseq->max_len, model->malloced_states, 2));
        //float back[4];
        for(i = 0; i < 4;i++){
                model->back[i] = 0.0f;
        }
        for(i = 0; i < numseq;i++){
                for(j = 0; j < iseq->len[i];j++){
                        iseq->labels[i][j] = (int)rk_interval(model->expected_K -1,   &model->rndstate)+2 ; // +2 for start and end state...
                        model->back[(int) iseq->seq[i][j]]++;
                }
        }
	
        /*for(i =0;i < iseq->num_seq;i++){
          for(j = 0; j < iseq->len[i];j++){
          fprintf(stdout,"  %c ","ACGT"[(int)iseq->seq[i][j]]);
          }
          fprintf(stdout,"\n");
          for(j = 0; j < iseq->len[i];j++){
          fprintf(stdout,"%3d ",iseq->labels[i][j]);
          }
          fprintf(stdout,"\n");
          }
        */
	
        
	
        for(i = 0; i < 4;i++){
                sum += model->back[i];// 0.0f;
        }
        for(i = 0; i < 4;i++){
                model->back[i] /= sum;
                model->back[i] = prob2scaledprob(model->back[i]);
        }
	
        //DPRINTF2("GOT background.");
        clear_counts(model);
        add_random_counts(model,100); 
        //fill_counts(iseq, model);
        //trim_counts(model);
        for(i = 0; i < 5;i++){
                RUN(iHmmHyperSample(model,10));
        }
        sample_counts(model);
	
        //print_emisson_matrix(model);
        //print_transistion_matrix(model);
        //print_hyper_parameters(model);
        //exit(0);
        float s0,s1,s2;
        for(iter = 0; iter <= model->numb + (model->nums -1)  *model->numi; iter++){//(numb + (nums-1)*numi);iter++){	
                for(i = 0; i < numseq;i++){
                        RUN(pgas_sample(model,pgas, iseq, i));
                }
                //	DPRINTF2("SAMPLING DONE.");
                clear_counts(model);
                remove_unused_states_smc(model, iseq);
                //add_random_counts(model,1); 
                fill_counts(iseq, model);
                //trim_counts(model);
		
                RUN(iHmmHyperSample(model,20));
                sample_counts(model);
                fprintf(stderr,"Iteration %4d: K = %3d, alpha0 = %f, gamma = %f\n",iter, model->K+1, model->alpha0, model->gamma);
		
                //model->alpha0 *= 0.99;//(float) total_len*1;
                //model->gamma *= 0.99;
                //print_prior(model);
                //print_transistion_matrix(model);
                //print_emisson_matrix(model);
                if (iter >= model->numb && ((iter-model->numb) % model->numi) == 0){
                        model->collect_alpha += model->alpha0;
                        model->collect_gamma += model->gamma;
                }
        }
	
        model->collect_gamma = model->collect_gamma / (float) model->nums;
        model->collect_alpha = model->collect_alpha / (float) model->nums;
	
        s0 = 0.0f;
        s1 = 0.0f;
        s2 = 0.0f;
	
        for(i = 0 ; i< iseq->num_seq;i++){
                sum = prob2scaledprob(1.0);
                for(j = 0; j < iseq->len[i];j++){
                        sum = sum + model->back[(int) iseq->seq[i][j]] + prob2scaledprob(1.0 /(float) iseq->len[i]);
                }
                sum += prob2scaledprob(1.0 - (1.0 / (float)iseq->len[i]));
                //fprintf(stdout,"%f %f \n",sum,iseq->score[i] );
                s0++;
                s1 += (iseq->score[i] - sum);
                s2 += (iseq->score[i] - sum)*  (iseq->score[i] - sum);
        }
        fprintf(stderr,"Log likelihood: %f %f	%f	%f\n",s1/s0, sqrt(  (s0 * s2 - pow(s1,2.0))   /  (  s0 *(s0-1.0) )),model->collect_alpha,model->collect_gamma);
        
        for(i = 0; i <= model->K;i++){
                model->sumM[i] = 0;
        }
        for(i = 0; i < iseq->num_seq;i++){
                for(j= 0; j <  iseq->len[i];j++){
                        model->sumM[iseq->labels[i][j]]++;
                }
        }
	
	
	
        //for(i = 0; i < numseq;i++){
        //	MFREE(sequences[i]);//, sizeof(char)* (len+1));
        //}
	
        //MFREE(sequences);//, numseq, sizeof(char*) );
        //free_model(model);
        print_emisson_matrix(model);
        
        free_ihmm_seq(iseq);
	
        free_pgas(pgas);
	
	
        return OK;
ERROR:
        return FAIL;
}



int pgas_sample(struct iHMM_model* model, struct pgas* pgas,struct ihmm_sequences* iseq,int num)
{
        int i,j,c,g;
        char* seq = NULL;
        int* lab = NULL;
        int len;
	
        int letter = -1;
        float tmp_x = 0.0;
        float* x = pgas->x;
        float* y = pgas->y;
        float** emission = model->emission;
        float** transition = model->transition;
        float sum = 0.0;
	
        float sum2 = 0.0;
	
        float marginal_likelihood = 0.0;
        float m_sum = 0.0;
	
        seq = iseq->seq[num];
        lab =  iseq->labels[num];
        len =  iseq->len[num];
	
	
        /*
          Calculate the probability of the previous path using the current
          transition and emission probabiloities. Note: this is done backward so
          we can concatenate sampled paths with the old path for weighting
          ancestors later.
        */
        pgas->previous_path_score[len] = prob2scaledprob(1.0);
	
        letter = (int)seq[len-1];
        pgas->previous_path_score[len-1] = prob2scaledprob(model->emission[lab[len-1] ][ letter]) + prob2scaledprob( model->transition[lab[len-1]][iHMM_STOP_STATE]  );
        for(j = len -2;j >= 0; j--){
                letter = (int)seq[j];
                pgas->previous_path_score[j] =
                        pgas->previous_path_score[j+1] + // to next state...
                        prob2scaledprob(model->transition[lab[j]][lab[j+1]]) + // transition
                        prob2scaledprob(model->emission[lab[j]][letter]); // emission
        }
	
        //pretty_print_vector_float("Previous path scores:",pgas->previous_path_score,len,1);
	
	
	
        // figure out probability of seeing y[0] in any state
	
        //sample from this ...
        sum = 0.0;
        letter = (int)seq[0];
        //P(xj | y1);
        x[iHMM_START_STATE] = 0.0;
        for(j = 1;j <= model->K;j++){
                /*if(j == 0){
                  x[j] = emission[j][letter];
                  }else{
                  x[j] = 0.0;
                  }*/
		
                x[j] = emission[j][letter]*   model->transition[iHMM_START_STATE][j];// model->prior[j];// (1.0f/ (float) (model->K+1)) ;// * 1.0f/ (float) model->K;
                sum += x[j];
                //sum = logsum(sum, x[j]);
        }
        for(j = 1;j <= model->K;j++){
                // q ( xj | y1)
                x[j] /=sum;
                //DPRINTF3("%f", x[j]);

                //y[j] = 0.0;
        }
	
	
        //fprintf(stdout,"letter %d: %d\n",0, letter);
	
        //pretty_print_vector_float("P (Xj | y1) - to be sampled from...",x,model->K+1,0);
	
        // sample particels ;
	
        marginal_likelihood = prob2scaledprob(1.0);
        m_sum = 0.0;
        sum = 0.0;
        for(j = 0; j < pgas->num_particles;j++){
                c = -1;
                RUN( select_random(x,model->K+1,&c,  rk_double(&model->rndstate)));
                //expand if we select infinityghost...
		
		
                pgas->particle_ancestry[0][j] = c;
                m_sum += x[c];
		
                //	DPRINTF3("%d	selected",c);
                /*if(c ==0){
                  tmp_x =emission[c][letter];
                  }else{
                  tmp_x = 0.0;
                  }*/
                tmp_x = emission[c][letter] *  model->transition[iHMM_START_STATE][c];//   model->prior[c];// * (1.0f/ (float) (model->K+1)) *
		
                //	DPRINTF3("prob: prior %f * emission %f = %f\n",  (1.0f/ (float)( model->K+1)) ,emission[c][letter], tmp_x);
                pgas->particle_path_prob[j] = prob2scaledprob(tmp_x);
		
                //weigth xj = sum prob of getting there for each particle...
                y[j] = tmp_x;
                sum += tmp_x;
                //	DPRINTF3("sum : %f \n",sum);
		
        }
        marginal_likelihood += prob2scaledprob(m_sum / pgas->num_particles);
        //What is the weigth of state 1,2 ... K ?
	
        for(j = 0;j  <pgas->num_particles;j++){
                y[j] =    prob2scaledprob(y[j] /  sum);
		
                ///weighth of state j according to particles is :
                //1 ) probability of getting there given y1
                // 2) divided by all particess
                //x[j] =(1.0f/ (float) (model->K+1)) *emission[j][(int)seq[0]] / sum ;
        }
	
	
        //pretty_print_vector_int("Initial particles:",pgas->particle_ancestry[0], pgas->num_particles);
        //pretty_print_vector_float("wj: prior (Xj) * emission (y1 | Xj ) / q( Xj| y1 )",y, pgas->num_particles,1);
        //pretty_print_vector_float("Path prob:",pgas->particle_path_prob, pgas->num_particles,1);
	
        //exit(0);
        for (i = 1; i <  len;i++){
                //get ancestor weights
                letter = (int)seq[i];
                //	fprintf(stdout,"letter %d: %d\n",i, letter);
                //for(j = 0;j <= model->K;j++){
                //	pgas->ancestor_weight[j] = 0.0;
                //}
                //sum = 0.0;
                for(j = 0;j < pgas->num_particles;j++){
			
                        c =pgas->particle_ancestry[i-1][j] & 0xFFFFFF; // current state of path [j]..
                        //if(pgas->ancestor_weight[c] == 0.0){
                        //DPRINTF3("State: %d\n",c);
                        tmp_x =
                                pgas->particle_path_prob[j] + // path up until here
                                prob2scaledprob(transition[c][ lab[i]] ) + // transiton to state in previous path
                                pgas->previous_path_score[i]; // score of old path from j until len
                        //DPRINTF3("%f %f %f : %f\n",pgas->particle_path_prob[j],prob2scaledprob(transition[c][ lab[i]] ),pgas->previous_path_score[i],tmp_x);
                        //DPRINTF3("%d: %e\n",c,scaledprob2prob(tmp_x - pgas->particle_path_prob[j]));
                        if(tmp_x > prob2scaledprob(0.0)){
                                pgas->ancestor_weight[j]  = y[j]  + (tmp_x-  pgas->particle_path_prob[j] ); // whole path divided by bit so far...
                        }else{
                                pgas->ancestor_weight[j] = prob2scaledprob(0.0);
                        }
                        //	sum +=pgas->ancestor_weight[c];
                        //}
                        //fprintf(stdout,"%d %f\n",j, scaledprob2prob(tmp_x -pgas->particle_path_prob[j]));
                }
                sum = prob2scaledprob(0.0);
                for(j = 0;j <pgas->num_particles;j++){
                        //	DPRINTF3("sum:%f next %f\n", sum,pgas->ancestor_weight[j]  );
                        sum = logsum(sum,  pgas->ancestor_weight[j]);
                        //sum +=pgas->ancestor_weight[j];
                }
		
                for(j = 0;j <pgas->num_particles;j++){
                        pgas->ancestor_weight[j] = pgas->ancestor_weight[j] -  sum;
                        pgas->ancestor_weight[j] = scaledprob2prob(pgas->ancestor_weight[j]);
                }
                //sum = 0.0;
                //for(j = 0;j <pgas->num_particles;j++){
                //	sum += pgas->ancestor_weight[j];
                //sum +=pgas->ancestor_weight[j];
                //}
                //for(j = 0;j <pgas->num_particles;j++){
                //	pgas->ancestor_weight[j] /= sum;
                //sum +=pgas->ancestor_weight[j];
                //}
		
                //ancestor_weight
                //	pretty_print_vector_float("ancestor weigth = w j *  (prob path new 0->i concat  old i > len ) / prob new path 0->i ",pgas->ancestor_weight, pgas->num_particles,0);
		
                //exit(0);
                // sample ancestors...
                sum = 0.0f;
                //for(j = 0;j < pgas->num_particles;j++){
                //	y[j] = 0.0;
                //	pgas->ancestors[j] = -1;
                //}
                //get ancestors
                for(j = 0;j < pgas->num_particles;j++){
                        c = -1;
                        RUN( select_random(pgas->ancestor_weight,pgas->num_particles,&c , rk_double(&model->rndstate)));
                        pgas->ancestors[j] = c; // ancestor j suggest to trace back to state c
			
                }
                // given than ancestor of partile j is state c what is the probability of ending up in x 1, x2 ,x3....
                //	pretty_print_vector_int("Picked ancestors..",pgas->ancestors, pgas->num_particles);
                sum2 = 0.0;
                m_sum = 0.0f;
                for(j = 0;j < pgas->num_particles;j++){
                        for(c = 0; c <= model->infinityghost;c++){
                                x[c] = 0;
                        }
                        sum = 0.0;
			
			
                        g = pgas->particle_ancestry[i-1][pgas->ancestors[j]] & 0xFFFFFF ;
                        //		DPRINTF3("Ancestor of %d is a %d state\n", j,g);
                        for(c = 2; c <= model->K ;c++){
                                x[c] = transition[g][c]  *emission[c][letter];
                                sum += x[c];
                        }
                        c = model->infinityghost;
                        x[c] =  rk_gamma(&model->rndstate, EMISSION_H, 1.0);
                        x[c]  = x[c]  / ( x[c]  + rk_gamma(&model->rndstate, EMISSION_H, 1.0)+ rk_gamma(&model->rndstate, EMISSION_H, 1.0)+ rk_gamma(&model->rndstate, EMISSION_H, 1.0));
                        x[c]= x[c] * transition[g][c];
                        sum += x[c];
			
                        for(c = 0; c <= model->infinityghost;c++){
                                x[c] /= sum;
                                //	fprintf(stdout,"%f ",x[c]);
                        }
                        //fprintf(stdout,"\n");
                        //fflush(stdout);
                        //		pretty_print_vector_float("particle weigth: ",x,model->K+1,0);
			
                        // which state are we going to from  ancestor ?
                        c = -1;
                        RUN( select_random(x,model->infinityghost+1,&c, rk_double(&model->rndstate) ));
                        //		fprintf(stdout,"Picked: %d\n",c);
                        m_sum = m_sum +x[c];
                        //c =model->infinityghost;
                        if(c == model->infinityghost){
				
                                //	fprintf(stderr,"Generated new STATE!!!!\n");//
                                RUN(add_state_to_model(model));
                                RUN(resize_pgas(pgas, model->malloced_states));
				
                                //print_transistion_matrix(model);
				
                                transition = model->transition;
                                emission = model->emission;
                                x = pgas->x;
                                y = pgas->y;
                                // resize pgas...
				
                        }
			
                        pgas->particle_ancestry[i][j] = (pgas->ancestors[j] << 24) | c;
			
                        g = pgas->particle_ancestry[i-1][pgas->ancestors[j]] & 0xFFFFFF;
                        //fprintf(stdout,"from:: %d	to %d %d\n",g,c,pgas->ancestors[j]);
                        //fflush(stdout);
                        //		fprintf(stdout,"Calc path prob\n");
			
                        //fprintf(stdout,"%f	- so far\n", scaledprob2prob(pgas->particle_path_prob[pgas->ancestors[j]]));
                        //fflush(stdout);
                        //fprintf(stdout,"%f	- transition\n",transition[g][c]);
                        //fflush(stdout);
                        //fprintf(stdout,"%f	- emission\n",emission[c][letter]);
                        //fflush(stdout);
			
                        pgas->tmp_particle_path_prob[j] =pgas->particle_path_prob[pgas->ancestors[j]] + prob2scaledprob(transition[g][c]  *emission[c][letter]);
			
			
                        //if(pgas->ancestors[j]!= -1){
                        g =pgas->particle_ancestry[i-1][pgas->ancestors[j]] & 0xFFFFFF;
                        //r(c = 0; c < pgas->num_particles;c++){.
                        // weight - transition from ancestor to state 'c' divided by prob of landing htere given sequence[i] ...
                        sum2 += x[c];
                        y[j] = transition[g][c]  *emission[c][letter]; //  / x[c];
                        //sum += x[j];
                        //}
                        //fprintf(stdout,"I is equal to %d\n",i);
                        //print_transistion_matrix(model);
			
                }
		
                //if(i == 8){
                //	exit(0);
                //}
                marginal_likelihood += prob2scaledprob(m_sum / pgas->num_particles)  ;
                sum = 0.0;
                for(j = 0;j < pgas->num_particles;j++){
                        pgas->particle_path_prob[j] = pgas->tmp_particle_path_prob[j];
                        //	fprintf(stdout,"%f ",pgas->particle_path_prob[j] );
                        y[j] /= sum2;
                        sum += y[j];
                }
                //fprintf(stdout,"\n");
                //fflush(stdout);
                for(j = 0;j < pgas->num_particles;j++){
                        y[j] /= sum;
                }
		
                //	pretty_print_vector_float("particle weigth: ",y, pgas->num_particles,0);
                //pretty_print_vector_float("Path prob:",pgas->particle_path_prob, pgas->num_particles,0);
                //exit(0);
		
        }
        // transition to stop !!!!
	
        i = len;
        for(j = 0;j < pgas->num_particles;j++){
		
                c =pgas->particle_ancestry[i-1][j] & 0xFFFFFF; // current state of path [j]..
                //if(pgas->ancestor_weight[c] == 0.0){
                //DPRINTF3("State: %d\n",c);
                tmp_x =
                        pgas->particle_path_prob[j] + // path up until here
                        prob2scaledprob(transition[c][ iHMM_STOP_STATE] ) + // transiton to state in previous path
                        pgas->previous_path_score[i]; // score of old path from j until len
                //DPRINTF3("%f %f %f : %f\n",pgas->particle_path_prob[j],prob2scaledprob(transition[c][ lab[i]] ),pgas->previous_path_score[i],tmp_x);
                //DPRINTF3("%d: %e\n",c,scaledprob2prob(tmp_x - pgas->particle_path_prob[j]));
                if(tmp_x > prob2scaledprob(0.0)){
                        pgas->ancestor_weight[j]  = y[j]  + (tmp_x-  pgas->particle_path_prob[j] ); // whole path divided by bit so far...
                }else{
                        pgas->ancestor_weight[j] = prob2scaledprob(0.0);
                }
                //	sum +=pgas->ancestor_weight[c];
                //}
                //fprintf(stdout,"%d %f\n",j, scaledprob2prob(tmp_x -pgas->particle_path_prob[j]));
        }
        sum = prob2scaledprob(0.0);
        for(j = 0;j <pgas->num_particles;j++){
                //	DPRINTF3("sum:%f next %f\n", sum,pgas->ancestor_weight[j]  );
                sum = logsum(sum,  pgas->ancestor_weight[j]);
                //sum +=pgas->ancestor_weight[j];
        }
	
        for(j = 0;j <pgas->num_particles;j++){
                pgas->ancestor_weight[j] = pgas->ancestor_weight[j] -  sum;
                pgas->ancestor_weight[j] = scaledprob2prob(pgas->ancestor_weight[j]);
        }
	
        for(j = 0;j < pgas->num_particles;j++){
                c = -1;
                RUN( select_random(pgas->ancestor_weight,pgas->num_particles,&c , rk_double(&model->rndstate)));
                pgas->ancestors[j] = c; // ancestor j suggest to trace back to state c
		
        }
	
        for(j = 0;j < pgas->num_particles;j++){
                g = pgas->particle_ancestry[i-1][pgas->ancestors[j]] & 0xFFFFFF;
                pgas->tmp_particle_path_prob[j] =pgas->particle_path_prob[pgas->ancestors[j]] + prob2scaledprob(transition[g][iHMM_STOP_STATE]);
        }
        for(j = 0;j < pgas->num_particles;j++){
                pgas->particle_path_prob[j] = pgas->tmp_particle_path_prob[j];
        }
	
        //fprintf(stdout,"Marginal: %f\n", marginal_likelihood );
	
	
        iseq->score[num] = marginal_likelihood;
        //backtrace
        c = -1;
        sum = prob2scaledprob(0.0);
	
        for(j = 0;j < pgas->num_particles;j++){
                sum = logsum(sum, pgas->particle_path_prob[j]);
                //if(pgas->particle_path_prob[j] >sum){
                //	sum =pgas->particle_path_prob[j];
                //	c = j;
                //}
        }
        sum2 = 0.0;
        for(j = 0;j < pgas->num_particles;j++){
                pgas->particle_path_prob[j] = scaledprob2prob(pgas->particle_path_prob[j] -sum);
                sum2 +=pgas->particle_path_prob[j];
        }
	
        for(j = 0;j < pgas->num_particles;j++){
                pgas->particle_path_prob[j] /= sum2;
        }
        c = -1;
        RUN( select_random(pgas->particle_path_prob, pgas->num_particles,&c, rk_double(&model->rndstate)));
	
	
        //fprintf(stdout,"SCORE:%f\n",sum);
	
        j = -1;
        sum = prob2scaledprob((1.0f/ (float) (model->K+1)));
        //fprintf(stdout,"Sum:%f\n",sum);
        for(i = len-1;i >= 0 ;i--){
                g =pgas->particle_ancestry[i][c]& 0xFFFFFF;
		
                iseq->labels[num][i] = g;
                //fprintf(stdout,"%d %d	\n%f\n",pgas->particle_ancestry[i][c]& 0xFFFFFF, pgas->particle_ancestry[i][c] >> 24,scaledprob2prob(sum));
                if(j == -1){
                        sum = sum + prob2scaledprob(emission[g][(int)seq[i]]);
                }else{
                        sum = sum + prob2scaledprob(emission[g][(int)seq[i]]  *transition[g][j] );
                }
                //fprintf(stdout,"%d Sum:%f\n",i, sum);
                j =g;
                c = pgas->particle_ancestry[i][c] >> 24;
        }
        return OK;
ERROR:
	
        return FAIL;
}

int remove_unused_states_smc(struct iHMM_model* model, struct ihmm_sequences* iseq)
{
        int i,j;
	
        float sum = 0.0;
        for(i = 0; i <= model->K;i++){
                model->sumM[i] = 0;
        }
        model->sumM[iHMM_START_STATE] = 10000;
        model->sumM[iHMM_STOP_STATE] = 10000;
	
        for(i = 0; i < iseq->num_seq;i++){
                for(j= 0; j <  iseq->len[i];j++){
                        model->sumM[iseq->labels[i][j]]++;
                }
        }
        j = 0;
        for(i = 0; i <= model->K;i++){
                if(model->sumM[i] != 0){
                        //if(dyn->state_used[i]  > (double) sb->numseq * 0.25){
			
                        model->beta[j] = model->beta[i];
                        model->sumN[i] = j;
                        //dyn->ID[i] = j;
                        j++;
                }else{
                        model->sumN[i] = j;
                        //dyn->ID[i] = j;
                        sum += model->beta[i];
                }
                //dyn->ID[i] = j;
                //fprintf(stderr,"%d %d %d\n",i,model->sumM[i],model->sumN[i] );
        }
        //fprintf(stderr,"%d	%f\n",j,sum);
        model->beta[j] = sum;
        for(i = j+1; i <= model->K;i++){
                model->beta[i] = 0;
        }
        model->K = j-1;
        model->infinityghost = j;
        for(i = 0; i < iseq->num_seq;i++){
                for(j= 0; j <  iseq->len[i];j++){
                        iseq->labels[i][j] = model->sumN[iseq->labels[i][j] ];
                        //fprintf(stderr,"Lab: %d  (seq:%d pos:%d \n",iseq->labels[i][j] ,i,j);
                        //sb->sequences[i]->labels[j] = dyn->ID[sb->sequences[i]->labels[j]];
                        //	sequence->labels[i] = id[sequence->labels[i]];
			
                        //		KSL_DASSERT3((sb->sequences[i]->labels[j] > 0  ));
                }
        }
        //RUN(resize_ (model, <#int num_states#>, <#int len#>)
        //resize_sample(sample);
        return OK;
}



struct pgas* init_pgas(int max_len, int malloced_states,int num_particles)
{
        struct pgas* pgas = NULL;
	
        MMALLOC(pgas, sizeof(struct pgas));
        pgas->particle_ancestry = NULL;
        pgas->particle_path_prob = NULL;
        pgas->tmp_particle_path_prob = NULL;
        pgas->previous_path_score = NULL;
        pgas->x = NULL;
        pgas->y = NULL;
        pgas->ancestors = NULL;
        pgas->ancestor_weight = NULL;
        pgas->num_particles = num_particles;
        pgas->buff_size = (malloced_states << 1) + num_particles;
        pgas->particle_ancestry = malloc_2d_int(pgas->particle_ancestry,max_len, pgas->num_particles , -1);
        MMALLOC(pgas->particle_path_prob, sizeof(float) * pgas->num_particles);
        MMALLOC(pgas->tmp_particle_path_prob, sizeof(float) * pgas->num_particles);
        MMALLOC(pgas->previous_path_score,sizeof(float) * (max_len+1));
        MMALLOC(pgas->x, sizeof(float) * pgas->buff_size);
        MMALLOC(pgas->y, sizeof(float) * pgas->buff_size);
        MMALLOC(pgas->ancestors, sizeof(int) * pgas->num_particles );
        MMALLOC(pgas->ancestor_weight, sizeof(float) * pgas->num_particles);
        return pgas;
ERROR:
        free_pgas(pgas);
        return NULL;
}

int resize_pgas(struct pgas* pgas, int malloced_states)
{
        int x = (malloced_states << 1)+pgas->num_particles;
        if(x > pgas->buff_size ){
                pgas->buff_size = x;
                MREALLOC(pgas->x, sizeof(float) * pgas->buff_size);
                MREALLOC(pgas->y, sizeof(float) * pgas->buff_size);
        }
        return OK;
ERROR:
        return FAIL;
}

void free_pgas(struct pgas* pgas)
{
        if(pgas){
                free_2d((void**) pgas->particle_ancestry);
                MFREE(pgas->particle_path_prob);
                MFREE(pgas->tmp_particle_path_prob);
                MFREE(pgas->previous_path_score);
                MFREE(pgas->x);
                MFREE(pgas->y);
                MFREE(pgas->ancestors);
                MFREE(pgas->ancestor_weight);
                MFREE(pgas);
        }
}

int iHmmSampleBeam(struct ihmm_sequences* iseq,struct iHMM_model* model)
{
        //rk_randomseed( &rndstate);
	
        struct dp* dp = NULL;
        int i,j,iter;
        float min_u = 1.0;
        float max_into_infinity = 0.0;
        float old_average = 0.0;
        float score = 0;
        RUN(start_iHMM_model(model, model->expected_K));
	
        RUNP(dp = init_dp(model->infinityghost, iseq->max_len ));
	
        RUN(clear_counts(model));
        RUN(fill_counts(iseq, model));
	
        RUN(print_trans_count_matrix(model));
	
        for(i = 0; i < 5;i++){
                //print_hyper_parameters(model);
                RUN(iHmmHyperSample(model,100));
        }
	
        RUN(sample_counts(model));
        print_transistion_matrix(model);
        print_emisson_matrix(model);
	
        for(iter = 0; iter <= model->numb + (model->nums -1)  *model->numi; iter++){//(numb + (nums-1)*numi);iter++){
                if( iter < model->numb){
                        //KSL_DPRINTF3(("new:%d old: %d\n", sample->K,sample->malloced_K ));
                        fprintf(stderr,"Iteration %d: K = %d, alpha0 = %f, gamma = %f\n",iter, model->K+1, model->alpha0, model->gamma);
                }
                RUN(clear_counts(model));
		
                // sequence_counter = 0;
                // clear_counts(sample);
                //
                RUN(set_u(model, iseq, &min_u));
		
                max_into_infinity = 0.0;
		
                for(i = 0; i <= model->K;i++){
                        if(model->transition[i][model->K] > max_into_infinity){
                                max_into_infinity = model->transition[i][model->K];
                        }
                }
                //fprintf(stdout,"%f minu : %f\n", max_into_infinity,min_u );
                while(max_into_infinity > min_u){
                        RUN(add_state_to_model(model));
                        //RUN(print_transistion_matrix(model),"print transtion matrix failed.");
                        max_into_infinity = 0.0;
			
                        for(i = 0; i <= model->K;i++){
                                if(model->transition[i][model->K] > max_into_infinity){
                                        max_into_infinity = model->transition[i][model->K];
                                }
                        }
                        //fprintf(stdout,"%f minu : %f\n", max_into_infinity,min_u );

                }
                RUN(resize_dp(dp, model->infinityghost, iseq->max_len));
                for(i = 0;i < model->infinityghost;i++){
                        dp->state_used[i] = 0;
                        dp->ID[i] = 0;
                }
		
                float s0,s1,s2;
                s0 = 0.0;
                s1 = 0.0;
                s2 = 0.0;
		
                //RUN(clear_counts(model),"clear counts failed.");
                for(i = 0; i < iseq->num_seq;i++){
                        RUN(dyn_prog(dp, model, iseq,i));
                        score = log(1.0);
                        //ASSERT(model->prior[iseq->labels[i][0]] > 0,"prior is zero");
                        score = score + log(model->prior[iseq->labels[i][0]]);
                        score = score + log(model->emission[iseq->labels[i][0]][(int)iseq->seq[i][0] ]);
                        for(j =1; j < iseq->len[i];j++){
                                ASSERT(model->transition[iseq->labels[i][j-1]][iseq->labels[i][j]] > 0,"transition zero %d - %d:%f u: %f %f \n",j-1,j, model->transition[iseq->labels[i][j-1]][iseq->labels[i][j]] , iseq->u[i][j] , iseq->u[i][j-1]);
                                ASSERT(model->emission[iseq->labels[i][j]][(int)iseq->seq[i][j] ]> 0,"emission zero\n");
                                score = score + log(model->transition[iseq->labels[i][j-1]][iseq->labels[i][j]]);
                                score = score + log(model->emission[iseq->labels[i][j]][(int)iseq->seq[i][j] ]);
                        }
			
                        s0++;
                        s1 += score;
                        s2 += score * score;
                }
		
                fprintf(stderr,"%f	%f\n",s1/s0, sqrt(  (s0 * s2 - pow(s1,2.0))   /  (  s0 *(s0-1.0) )));
		
                RUN(remove_unused_states(dp, model,iseq));
		
                RUN(fill_counts(iseq, model));
	
                RUN(iHmmHyperSample(model,100));
		
                RUN(sample_counts(model));
		
                if(iter % 10 == 0){
                        if(fabs(old_average - (s1/s0)) < 0.01){
                                //	break;
                        }
                        old_average =s1/s0;
                }
                if (iter >= model->numb && ((iter-model->numb) % model->numi) == 0){
                        model->collect_alpha += model->alpha0;
                        model->collect_gamma += model->gamma;
                }

        }
        model->collect_gamma = model->collect_gamma / (float) model->nums;
        model->collect_alpha = model->collect_alpha / (float) model->nums;
        RUN(free_dp(dp));
        return OK;
	
ERROR:
        return FAIL;
}

int dyn_prog(struct dp* dyn,struct iHMM_model* model, struct ihmm_sequences* iseq,int num)
{
	
        float** dp = NULL;
	
        float* u = NULL;
	
        int i,j,c,old_j;
	
	
	
        float sum;
        int* state_used = NULL;
        int* id = NULL;
        int*label = NULL;
        RUN(clear_dp(dyn));
	
	
	
        char* seq = NULL;
        int len =  iseq->len[num];
        //sequence->len;
	
        state_used = dyn->state_used;
        id = dyn->ID;
	
        dp = dyn->dp_matrix;
	
	
        label = iseq->labels[num];
        u =  iseq->u[num];// sequence->u;
        seq = iseq->seq[num];
	
        //seq = sequence->seq;
	
        for(i = 0; i <= model->K; i++){
                if(model->prior[i] > u[0]){
                        //if(model->transition[0][i] >u[0]){
                        dp[0][i] = 1.0;
                }
        }
        sum = 0.0;
        for(i = 0; i <= model->K; i++){
                dp[0][i] = dp[0][i] * model->emission[i][(int)seq[0]];
                sum +=dp[0][i] ;
        }
        for(i = 0; i <= model->K; i++){
                dp[0][i] = dp[0][i] / sum;
        }
        for(i = 1; i < len;i++){
                sum = 0.0;
                for(j = 0; j <= model->K; j++){
                        for(c = 0; c <= model->K; c++){
                                //	if(model->transition[j][c] > u[i]){
                                dp[i][c] += dp[i-1][j] * (model->transition[j][c] > u[i]);
                                //}
                        }
                }
                for(j = 0; j  <= model->K;  j++){
                        dp[i][j] = dp[i][j] * model->emission[j][(int)seq[i]];
                        sum +=dp[i][j];
                }
                for(j = 0; j  <= model->K; j++){
                        dp[i][j] = dp[i][j] / sum;
                }
        }
        sum = 0.0;
        i = len-1;
        for(j = 0; j  <= model->K; j++){
                sum += dp[i][j];
        }
        if(sum != 0.0){
                i = len-1;
                sum = 0.0;
                for(j = 0; j  <= model->K; j++){
                        dp[i][j] = dp[i][j] + dp[i][j-1];
                }
                sum = dp[i][model->K];
                for(j = 0; j  <= model->K; j++){
                        dp[i][j] /= sum;
                }
                sum = rk_double(&model->rndstate);
                for(j = 0; j  <= model->K; j++){
                        if(sum <= dp[i][j]){
                                label[i] = j;
                                state_used[j]++;
                                break;
                        }
                }
		
                for(i= len-2;i >= 0;i--){
                        sum = 0.0;
                        old_j = 0;
                        for(j = 0; j <= model->K; j++){
                                if(model->transition[j][label[i+1]] > u[i+1]){
                                        dp[i][j] = dp[i][j] + sum;//dp[i][old_j];
                                        sum =dp[i][j] ;
                                        //old_j = j;
                                }
                        }
                        for(j = 0; j <= model->K; j++){
                                if(model->transition[j][label[i+1]] > u[i+1]){
                                        dp[i][j] /= sum;//dp[i][old_j] ;
                                }
                        }
			
                        sum = rk_double(&model->rndstate);
                        for(j = 0; j<= model->K; j++){
                                if(model->transition[j][label[i+1]] > u[i+1]){
                                        if(sum <= dp[i][j]){
                                                label[i] = j;
                                                state_used[j]++;
                                                break;
                                        }
                                }
                        }
                }
        }else{
                fprintf(stderr,"No path found\n");
        }
	
        return OK;
ERROR:
        return FAIL;
}

int remove_unused_states(struct dp* dyn,struct iHMM_model* model, struct ihmm_sequences* iseq)
{
        int i,j;
        j = 0;
        float sum = 0.0;
        for(i = 0; i <= model->K;i++){
                if(dyn->state_used[i] != 0){
                        //if(dyn->state_used[i]  > (double) sb->numseq * 0.25){
			
                        model->beta[j] = model->beta[i];
                        dyn->ID[i] = j;
                        j++;
                }else{
                        dyn->ID[i] = j;
                        sum += model->beta[i];
                }
                //dyn->ID[i] = j;
                //fprintf(stderr,"%d %d %d\n",i,dyn->state_used[i],j);
        }
	
        model->beta[j+1] = sum;
        for(i = j+2; i <= model->K;i++){
                model->beta[i] = 0;
        }
        model->K = j;
        model->infinityghost = j+1;
        for(i = 0; i < iseq->num_seq;i++){
                for(j= 0; j <  iseq->len[i];j++){
                        iseq->labels[i][j] = dyn->ID[iseq->labels[i][j] ];
                        //		fprintf(stderr,"Lab: %d  (seq:%d pos:%d \n",sb->sequences[i]->labels[j],i,j);
                        //sb->sequences[i]->labels[j] = dyn->ID[sb->sequences[i]->labels[j]];
                        //	sequence->labels[i] = id[sequence->labels[i]];
			
                        //		KSL_DASSERT3((sb->sequences[i]->labels[j] > 0  ));
                }
        }
        //RUN(resize_ (model, <#int num_states#>, <#int len#>)
        //resize_sample(sample);
        return OK;
}


int set_u(struct iHMM_model* model, struct ihmm_sequences* iseq, float* min_u)
{
        int i,j;
        float** u = iseq->u ;
        int** labels  = iseq->labels;
	
        *min_u  = 1.0;
	
        for(i = 0; i < iseq->num_seq;i++){
                ASSERT(model->prior[labels[i][0]] > 0.0  ,"weird: %f",model->prior[labels[i][0]]);
                u[i][0] = rk_double(&model->rndstate) * model->prior[labels[i][0]];

                //u[i][0] = rk_double(&rndstate) * model->transition[0][labels[i][0]];
                if(u[i][0] < *min_u){
                        *min_u = u[i][0];
                }
                for(j = 1; j < iseq->len[i];j++){
                        ASSERT(model->transition[labels[i][j-1]][labels[i][j]] > 0.0,"Weird   %d -> %d = %d.",labels[i][j-1],labels[i][j], model->transition[labels[i][j-1]][labels[i][j]]   );
                        u[i][j] = rk_double(&model->rndstate) * model->transition[labels[i][j-1]][labels[i][j]];
                        //fprintf(stderr,"%f %f   ( %d - > %d )\n", sample->transition[labels[i-1]][labels[i]], u[i],labels[i-1],labels[i] );
                        if(u[i][j] < *min_u){
                                *min_u = u[i][j];
                        }
                }
        }

        return OK;
ERROR:
        return FAIL;
}





int iHmmHyperSample(struct iHMM_model* model, int iterations)
{
        int i,j,c;
        float ialpha0 = model->alpha0;
        float igamma = model->gamma;
	
        int** N = model->trans_counts;
	
        int** M = model->M;
        int* sumM = model->sumM;
        int* sumN = model->sumN;
        int totalM = 0;
	
        float sum = 0.0;
	
        //MMALLOC(M, sizeof(int*) * (sample->K+2) );
        //MMALLOC(sumM, sizeof(int) * (sample->K+2) );
        //MMALLOC(sumN, sizeof(int) * (sample->K+2) );
        for(i = 0; i <= model->K;i++){
                //M[i] = NULL;
                sumM[i] = 0;
                sumN[i] = 0;
                //MMALLOC(M[i], sizeof(int) * (sample->K+2));
                for(j = 0; j  <= model->K;j++){
                        M[i][j] = 0;
                }
        }
	
        N = model->trans_counts;
	
        for(i = 1; i <= model->K;i++){
                for(j = 0; j <= model->K;j++){
                        if(N[i][j] == 0){
                                M[i][j] = 0;
                        }else{
                                for(c = 1;c <= N[i][j];c++){
                                        if(rk_double( &model->rndstate) < (ialpha0 *model->beta[j]) / (ialpha0 * model->beta[j] + (float)c -1.0)){
                                                M[i][j] = M[i][j] + 1; // number of times state i generated color j...
                                                totalM = totalM +1;
                                        }
                                }
                        }
                }
        }
	
	
        for(i = 1; i <= model->K;i++){
                for(j = 0; j <= model->K;j++){
                        sumM[j] += M[i][j];
                        sumN[i] += N[i][j];
                }
        }
        sum = 0.0;
        model->beta[0] = 0;
        for(i = 1; i <=  model->K;i++){
                model->beta[i] =rk_gamma(&model->rndstate, (float)sumM[i], 1.0);
                sum += model->beta[i];
        }
	
        model->beta[model->infinityghost] =  rk_gamma(&model->rndstate, igamma, 1.0);
        //fprintf(stdout,"BETA inf:%f\n",model->beta[model->infinityghost]  );
        sum += model->beta[model->infinityghost] ;
        for(i = 0; i <= model->infinityghost;i++){
		
                model->beta[i] /= sum;
                //	fprintf(stdout,"%f\n", model->beta[i] );
        }
	
        int iter;
        float* w = model->w;
        float* p = model->p;
        float* s = model->s;
        float sum_s = 0;
        float sum_w = 0;
	
	
        if(model->soft){
	 
                for(iter = 0;iter < iterations;iter++){
                        sum_s = 0.0;
                        sum_w = 0.0;
                        for(i = 0; i <= model->K;i++){
                                w[i]  = rk_beta(&model->rndstate,ialpha0+1.0,(float)sumN[i]);
                                p[i] = (float)sumN[i] / ialpha0;
                                p[i] = p[i] / (p[i]+1.0);
                                s[i]= rk_binomial(&model->rndstate, 1, p[i]);
                                sum_s +=s[i];
                                sum_w += log(w[i]);
                        }
                        ialpha0 = rk_gamma(&model->rndstate, model->alpha0_a + (float) totalM - sum_s, 1.0 / (model->alpha0_b -sum_w));
                }
                model->alpha0 = ialpha0;
        }
	
        float mu = 0.0;
        float pi_mu;
	
        if(model->soft){
                for(iter = 0;iter < iterations;iter++){
                        mu =  rk_beta(&model->rndstate, igamma+1.0, (float)totalM);
                        pi_mu = 1.0 / (1.0 + ((float)totalM * ( model->gamma_b - log(mu) )) / (model->gamma_a + (float)(model->infinityghost  ) -1.0));
                        if(rk_double( &model->rndstate) < pi_mu){
                                igamma = rk_gamma(&model->rndstate,model->gamma_a + (float)(model->infinityghost), 1.0 / (model->gamma_b - log(mu)));
                        }else{
                                igamma = rk_gamma(&model->rndstate, model->gamma_a + (float)(model->infinityghost)-1.0, 1.0 / (model->gamma_b - log(mu)));
                        }
                }
                model->gamma = igamma;
        }
	
	
        return OK;
}




int print_hyper_parameters(struct iHMM_model* model)
{
        int i;
        float sum = 0.0;
        fprintf(stdout,"BETA:\n");
        for(i = 0; i <= model->infinityghost;i++){
                sum+=model->beta[i];
                fprintf(stdout,"%f ",model->beta[i]);
        }
	
        fprintf(stdout,"\nSUM:%f\n",sum);
        fprintf(stdout,"ALPHA: %f\n",model->alpha0);
        fprintf(stdout,"GAMMA: %f\n",model->gamma);
        fflush(stdout);
        return OK;
}


int print_transistion_matrix(struct iHMM_model* model)
{
        int i,j;
        fprintf(stdout,"TRANSITION:\n");
        float sum = 0.0;
        for(i = 0; i <= model->infinityghost;i++){
                sum = 0.0;
                fprintf(stdout,"%d: ",i);
                for(j = 0;j <= model->infinityghost;j++){
                        sum +=model->transition[i][j];
                        fprintf(stdout,"%0.4f ",(model->transition[i][j]) );
                }
                fprintf(stdout,"	sum:%f\n",sum);
        }
        fprintf(stdout,"\n");
        fflush(stdout);
        return OK;
}

int print_trans_count_matrix(struct iHMM_model* model)
{
        int i,j;
        fprintf(stdout,"TRANSITION:\n");
        float sum = 0.0;
        for(i = 0; i <= model->infinityghost;i++){
                sum = 0.0;
                fprintf(stdout,"%d: ",i);
                for(j = 0;j <= model->infinityghost;j++){
                        sum +=model->transition[i][j];
                        fprintf(stdout,"%d ",(model->trans_counts[i][j]) );
                }
                fprintf(stdout,"	sum:%f\n",sum);
        }
        fprintf(stdout,"\n");
        fflush(stdout);
        return OK;
}



int print_prior(struct iHMM_model* model)
{
        int i;
        fprintf(stdout,"PRIOR:\n");
        float sum = 0.0;
        for(i = 0; i <= model->infinityghost;i++){
                sum+=model->prior[i];
                fprintf(stdout,"%0.4f ",(model->prior[i]) );
		
        }
        fprintf(stdout,"	sum:%f\n",sum);
        fflush(stdout);
        return OK;
}

int print_emisson_matrix(struct iHMM_model* model)
{
        int i,j;
        fprintf(stdout,"EMISSION:\n");
        float sum = 0.0;
        for(i = 0; i <= model->infinityghost;i++){
                sum = 0.0;
                fprintf(stdout,"%d: (used:%d)",i,model->sumM[i]);
                for(j = 0;j < model->L;j++){
                        sum +=model->emission[i][j];
                        fprintf(stdout,"%0.4f (%d)",(model->emission[i][j]),model->emit_counts[i][j] );
                }
                fprintf(stdout,"	sum:%f\n",sum);
        }
        fprintf(stdout,"\n");
        fflush(stdout);
        return OK;
}

int add_random_counts(struct iHMM_model* model,int count)
{
        int i,j;
	
        //for(i = 0; i <= model->K;i++){
        //	model->prior_counts[i] = 0.0;
        //}
	
        for(i = 0; i <= model->K;i++){
                for(j = 0; j <= model->K;j++){
                        model->trans_counts[i][j] += count;
                }
        }
        for(i = 0; i <= model->K;i++){
                for(j = 0; j < model->L;j++){
                        model->emit_counts[i][j] += count;
                }
        }
        return OK;
}


int fill_counts(struct ihmm_sequences* iseq,struct iHMM_model* model)
{
        int i,j;
	
        int** p = NULL;
	
        int** labels = iseq->labels;
	
        // get prior
        //for(i = 0; i < iseq->num_seq;i++){
        //	for(j =0; j < iseq->len[i];j++){
        //		model->prior_counts[labels[i][j]] += 1;
        //	}
        //}
	
        //get transistions...
        p = model->trans_counts;
        for(i = 0; i < iseq->num_seq;i++){
                p[iHMM_START_STATE][labels[i][0]] += 1;
                //p[0][labels[i][0]] += 1;
                for(j =1; j < iseq->len[i];j++){
                        p[labels[i][j-1]][labels[i][j]]++;
                }
		
                p[labels[i][iseq->len[i]-1]][iHMM_STOP_STATE] += 1;
        }
	
        //get emissions...
        p = model->emit_counts;
        for(i = 0; i < iseq->num_seq;i++){
                for(j =0; j < iseq->len[i];j++){
                        p[labels[i][j]][(int)iseq->seq[i][j]]++;
                }
		
        }
	
        return OK;
}

int trim_counts(struct iHMM_model* model)
{
        int i,j;
	
        for(i = 0;i <=model->K;i++){
                for(j = 0;j <=model->K;j++){
                        if(j <= i){
                                model->trans_counts[i][j] = 0;
                        }
                }
        }
        return OK;
}

int sample_counts(struct iHMM_model* model)
{
        int i,j;
        float sum = 0.0;
	
        //priot////
        /*for(j = 0; j <= model->K;j++){
          model->prior[j] = (float)model->prior_counts[j] + model->beta[j] * model->alpha0;
          model->prior[j]  = rk_gamma(&model->rndstate, model->prior[j] , 1.0);
          sum += model->prior[j] ;
          }
          model->prior[model->infinityghost] =  model->beta[model->infinityghost] * model->alpha0;
          model->prior[model->infinityghost]  = rk_gamma(&model->rndstate,model->prior[model->infinityghost], 1.0);
          sum += model->prior[j] ;
	
          for(j = 0;j <= model->infinityghost;j++){
          model->prior[j]  /= sum;
          }*/
        //START STATE !!!!
	
        i = iHMM_START_STATE;
        sum = 0.0;
        model->transition[i][iHMM_START_STATE] = 0.0f;
        model->transition[i][iHMM_STOP_STATE] = 0.0f;
	
        for(j = 2; j <= model->K;j++){
                model->transition[i][j] = (float)model->trans_counts[i][j] + model->beta[j] * model->alpha0;
                model->transition[i][j]  = rk_gamma(&model->rndstate, model->transition[i][j], 1.0);
                sum += model->transition[i][j] ;
        }
        model->transition[i][model->infinityghost] =  model->beta[model->infinityghost] * model->alpha0;
        model->transition[i][model->infinityghost]  = rk_gamma(&model->rndstate, model->transition[i][model->infinityghost], 1.0);
        sum += model->transition[i][model->infinityghost] ;
	
	
        for(j = 2;j <= model->infinityghost;j++){
                model->transition[i][j]  /= sum;
        }
        //STOP STATE
	
        i = iHMM_STOP_STATE;
        for(j = 0;j <= model->infinityghost;j++){
                model->transition[i][j] = 0.0f;
        }
	
	
	
        for(i = 2; i <= model->K;i++){
                sum = 0.0;
                //for(j = 0; j <=i;j++){
                //	model->transition[i][j]  =  0.0;//model->beta[j] * model->alpha0;
                //model->transition[i][j]  = rk_gamma(&model->rndstate, model->transition[i][j], 1.0);
                //	sum += model->transition[i][j] ;
                //}
                for(j = 1; j <= model->K;j++){
                        model->transition[i][j] = (float)model->trans_counts[i][j] + model->beta[j] * model->alpha0;
                        model->transition[i][j]  = rk_gamma(&model->rndstate, model->transition[i][j], 1.0);
                        sum += model->transition[i][j] ;
                }
                model->transition[i][model->infinityghost] =  model->beta[model->infinityghost] * model->alpha0;
                model->transition[i][model->infinityghost]  = rk_gamma(&model->rndstate, model->transition[i][model->infinityghost], 1.0);
                //if(i == model->K){
                //	model->transition[i][model->infinityghost]  =  1.0;
                //}
                //fprintf(stdout,"%f	%e\n",model->transition[i][model->infinityghost],rk_gamma(&model->rndstate, model->transition[i][model->infinityghost], 1.0) );
                sum += model->transition[i][model->infinityghost] ;

		
                for(j = 0;j <= model->infinityghost;j++){
                        model->transition[i][j]  /= sum;
                }
        }
	
        for(i = 2; i <= model->K;i++){
                sum = 0.0;
                for(j = 0; j < model->L;j++){
                        model->emission[i][j]  = model->emit_counts[i][j]  + EMISSION_H;
                        model->emission[i][j]  = rk_gamma(&model->rndstate, model->emission[i][j], 1.0);
                        sum +=model->emission[i][j];
                }
                for(j = 0;j <model->L;j++){
                        model->emission[i][j] /= sum;
                }
        }
        return OK;
}

int clear_counts(struct iHMM_model* model)
{
        int i,j;
	
        //for(i = 0; i <= model->K;i++){
        //	model->prior_counts[i] = 0.0;
        //}
	
        for(i = 0; i <= model->K;i++){
                for(j = 0; j <= model->K;j++){
                        model->trans_counts[i][j] = 0;
                }
        }
        for(i = 0; i <= model->K;i++){
                for(j = 0; j < model->L;j++){
                        model->emit_counts[i][j] = 0;
                }
        }
        return OK;
}



struct iHMM_model* init_iHMM_model(void)
{
        struct iHMM_model* model = NULL;
        MMALLOC(model, sizeof(struct iHMM_model ));
	
        model->L = 4;
        model->soft = 1;
        model->collect_alpha = 0.0;
        model->collect_gamma = 0.0;
        model->alpha0 = 0.0;
        model->gamma = 0.0;
        model->alpha0_a = 4.0;
        model->alpha0_b = 1.0;
        model->gamma_a =  3.0;
        model->gamma_b = 6.0;
        model->expected_K = 20;
        model->numb = 100;
        model->nums = 1;
        model->numi = 1;
        model->beta = NULL;     
	
        model->prior_counts = NULL;
        model->trans_counts = NULL;
        model->emit_counts = NULL;
	
        model->back = NULL;
        model->transition=NULL;
        model->emission=NULL;
        model->prior = NULL;
	
        model->w=NULL;
        model->p=NULL;
        model->s=NULL;
	
        model->M=NULL;
        model->sumM=NULL;
        model->sumN=NULL;
        model->tmp = NULL;
	
	
        rk_randomseed(&model->rndstate);

        return model;
ERROR:
        return NULL;
}


int add_state_to_model(struct iHMM_model* model)
{
	
        int tmp_old_infinityghost = 0;
        int tmp_old_k = 0;
        double sum = 0.0;
        float be;
        float bg;
        float a,b,pg,pe;
        int i,j;
	
        tmp_old_infinityghost = model->infinityghost;
        tmp_old_k = model->K;
	
	
	
	
	

	
        if(model->infinityghost+1 == model->malloced_states){
                model->malloced_states = model->malloced_states+10;
//		fprintf(stdout,"RESIZING: %d \n",model->malloced_states);
		
                MREALLOC(model->beta, sizeof(float) *model->malloced_states);
//		fprintf(stdout,"RESIZING: %d \n",model->malloced_states);
                MREALLOC(model->w, sizeof(float) *model->malloced_states);
                MREALLOC(model->p, sizeof(float) *model->malloced_states);
                MREALLOC(model->s, sizeof(float) *model->malloced_states);
		
                MREALLOC(model->prior, sizeof(float) *model->malloced_states);
                MREALLOC(model->prior_counts, sizeof(int) *model->malloced_states);
                MREALLOC(model->sumM, model->malloced_states*sizeof(int));
                MREALLOC(model->sumN, model->malloced_states*sizeof(int));
                MREALLOC(model->tmp,sizeof(double) *model->malloced_states);
		
		
		
                model->M = malloc_2d_int(model->M, model->malloced_states, model->malloced_states, 0);
                model->trans_counts = malloc_2d_int(model->trans_counts, model->malloced_states, model->malloced_states, 0);
                model->transition = malloc_2d_float(model->transition,model->malloced_states, model->malloced_states, 0.0);
		
                model->emit_counts = malloc_2d_int(model->emit_counts, model->malloced_states, model->L, 0);
                model->emission = malloc_2d_float(model->emission,model->malloced_states, model->L, 0.0);
        }
        model->K++;// num_states-1; // first state is 1 not zero ....
	
        model->infinityghost++;
        //fill new column in beta vector...
		
		
		
        //fill new row in transition....
        sum = 0.0;
	
        for(i = 1; i <= tmp_old_infinityghost;i++){
                model->tmp[i] = model->beta[i] * model->alpha0;
                model->tmp[i] = rk_gamma(&model->rndstate, model->tmp[i] , 1.0);
                sum+= model->tmp[i];
                //model->transition[tmp_old_infinityghost][i] = model->beta[i] * model->alpha0;
                //model->transition[tmp_old_infinityghost][i] = rk_gamma(&rndstate, model->transition[tmp_old_infinityghost][i], 1.0);
                //sum += model->transition[tmp_old_infinityghost][i];
        }
	
        ASSERT(sum != 0.0,"Weird new row has no data:%f %d\n",sum,tmp_old_infinityghost);
        model->transition[tmp_old_infinityghost][iHMM_START_STATE] = 0.0;
        for(i = 1; i <= tmp_old_infinityghost;i++){
                model->transition[tmp_old_infinityghost][i] =(float)( model->tmp[i]   / sum);
        }
	
        //fill new row in emission....
	
        sum = 0;
        for(i = 0; i < model->L;i++){
                model->emission[tmp_old_infinityghost][i] = rk_gamma(&model->rndstate, EMISSION_H, 1.0);
                sum +=model->emission[tmp_old_infinityghost][i];
        }
        for(i = 0; i < model->L;i++){
                model->emission[tmp_old_infinityghost][i] = model->emission[tmp_old_infinityghost][i] / sum;
        }
	
        //split last column of transistion matrix.... i.e.allow existing states to connect to new state...
	
        //first get beta for new column
	
        be = model->beta[tmp_old_infinityghost];
        bg = rk_beta(&model->rndstate, 1.0,model->gamma );
	
        model->beta[tmp_old_infinityghost] = bg*be;
        model->beta[tmp_old_infinityghost+1] = (1.0 - bg) *be;
	
        //now split prob in last columns...
	
        a = model->alpha0 * model->beta[tmp_old_infinityghost];
        b = 0.0;
        for(j = 0; j <= tmp_old_infinityghost;j++){
                b+=model->beta[j];
        }
        //fprintf(stdout,"	alpha:%f	a:%f b:%f\n",model->alpha0,a,b);
        b = model->alpha0 * (1.0 - b);
	
        for(i =0 ; i <= tmp_old_infinityghost;i++){
                if(a < 1e-2 || b < 1e-2){     // % This is an approximation when a or b are really small.
                        pg = rk_binomial(&model->rndstate, 1.0, a / (a+b));
                }else{
                        pg = rk_beta(&model->rndstate, a, b);
                }
                //f/printf(stdout,"PG:%f	%f	a:%f b:%f\n",pg,model->alpha0,a,b);

                pe = model->transition[i][tmp_old_infinityghost];
                model->transition[i][tmp_old_infinityghost] = pg * pe;
                model->transition[i][tmp_old_infinityghost+1] = (1.0-pg) * pe;
        }
	
        //prior
        if(a < 1e-2 || b < 1e-2){     // % This is an approximation when a or b are really small.
                pg = rk_binomial(&model->rndstate, 1.0, a / (a+b));
        }else{
                pg = rk_beta(&model->rndstate, a, b);
        }
	
        pe = model->prior[tmp_old_infinityghost];
        model->prior[tmp_old_infinityghost] = pg * pe;
        model->prior[tmp_old_infinityghost+1] = (1.0-pg) * pe;
	
        return OK;
ERROR:
	
        return FAIL;
}


int start_iHMM_model(struct iHMM_model* model, int num_states)
{
	
        int i;
	
	
        model->L = 4;
        model->K = num_states+1; // first state is 1 not zero ....
	
        model->infinityghost = num_states+2;
	
	
        model->malloced_states = num_states+3;
	
	
	
	
	
        //sample->S = NULL;
        //model->alpha0 = 0.0;
        //model->gamma = 0.0;
        model->beta = NULL;
        model->transition = NULL;
        model->emission = NULL;
        model->emit_counts = NULL;
        model->trans_counts = NULL;
	
        model->M = NULL;
        model->sumM = NULL;
        model->sumN = NULL;
	
        model->w = NULL;
        model->p = NULL;
        model->s = NULL;
	
	
        model->back = NULL;
	
        MCALLOC(model->back, model->L, sizeof(float));
        MMALLOC(model->beta, sizeof(float) *model->malloced_states);
        MMALLOC(model->w, sizeof(float) *model->malloced_states);
        MMALLOC(model->p, sizeof(float) *model->malloced_states);
        MMALLOC(model->s, sizeof(float) *model->malloced_states);
	
        MMALLOC(model->prior, sizeof(float) *model->malloced_states);
        MMALLOC(model->prior_counts, sizeof(int) *model->malloced_states);
        MMALLOC(model->tmp,sizeof(double) *model->malloced_states);
        MCALLOC(model->sumM, model->malloced_states,sizeof(int));
        MCALLOC(model->sumN, model->malloced_states,sizeof(int));
	
        model->M = malloc_2d_int(model->M, model->malloced_states, model->malloced_states, 0);
        model->trans_counts = malloc_2d_int(model->trans_counts, model->malloced_states, model->malloced_states, 0);
        model->transition = malloc_2d_float(model->transition,model->malloced_states, model->malloced_states, 0.0);
	
        model->emit_counts = malloc_2d_int(model->emit_counts, model->malloced_states, model->L, 0);
        model->emission = malloc_2d_float(model->emission,model->malloced_states, model->L, 0.0);
	
        if(model->soft){
                model->alpha0 = rk_gamma(&model->rndstate, model->alpha0_a,1.0 / model->alpha0_b);
                model->gamma = rk_gamma(&model->rndstate,model->gamma_a,1.0 / model->gamma_b);
        }
        for(i = 0; i <= model->infinityghost;i++){
                model->beta[i] = 1.0 / (float)(model->infinityghost+1);
        }
	
        return OK;
ERROR:
        free_iHMM_model(model);
        return FAIL;
}


int free_iHMM_model(struct iHMM_model* model)
{
        if(model){
                free_2d((void**) model->M);
                free_2d((void**) model->trans_counts);
                free_2d((void**) model->transition);
                free_2d((void**) model->emit_counts);
                free_2d((void**) model->emission);
		

                MFREE(model->back);
                MFREE(model->beta);
                MFREE(model->w);//=NULL;
                MFREE(model->p);//=NULL;
                MFREE(model->s);//=NULL;
	
                MFREE(model->tmp);
                MFREE(model->prior);
                MFREE(model->prior_counts);
	
                MFREE(model->sumM);//=NULL;
                MFREE(model->sumN);//=NULL;
                MFREE(model);
        }
        return OK;
}

int set_labels_based_on_alignment(struct ihmm_sequences* iseq,char** aln,int aln_len,int num_seq_in_aln)
{
        int i,j,c;
        for(i = 0; i < iseq->num_seq;i++){
                c = 0;
                for(j = 0; j < aln_len;j++){
                        //fprintf(stdout,"%c",aln[i][j]);
                        switch (aln[i][j]){
                        case 'A':
                        case 'C':
                        case 'G':
                        case 'T':
                                iseq->labels[i][c] = j;
                                c++;
                                break;
                        default:
                                break;
                        }
                }
        }
        return OK;
}


struct ihmm_sequences* init_ihmm_seq(char** sequences,int numseq)
{
        int i,j;
        struct ihmm_sequences* iseq = NULL;

	
        ASSERT(numseq !=0,"no sequences provided");
        MMALLOC(iseq, sizeof(struct ihmm_sequences));
        iseq->max_len = 0;
        iseq->len = NULL;
        iseq->u = NULL;
        iseq->seq = NULL;
        iseq->labels = NULL;
        iseq->score = NULL;
        iseq->num_seq = numseq;
	
        MMALLOC(iseq->len, sizeof(int) * numseq);
	
        for(i = 0;i < numseq;i++){
                if((j = (int)strlen(sequences[i])) > iseq->max_len){
                        iseq->max_len = j;
                }
                iseq->len[i] = j;
        }
        ASSERT(iseq->max_len > 1,"something wrong with sequence lengths");
	
        MMALLOC(iseq->score, sizeof(float) * numseq);
        iseq->labels = malloc_2d_int(iseq->labels , numseq, iseq->max_len, 0);
        iseq->u = malloc_2d_float(iseq->u, numseq,  iseq->max_len, 0);
        iseq->seq = malloc_2d_char(iseq->seq, numseq,  iseq->max_len, 0);
        for(i = 0;i < numseq;i++){
		
                iseq->score[i] = prob2scaledprob(0.0f);
                for(j = 0;j < iseq->len[i];j++){
                        switch(sequences[i][j]){
                        case 'A':
                        case 'a':
                                iseq->seq[i][j] = 0;
                                break;
                        case 'C':
                        case 'c':
                                iseq->seq[i][j] = 1;
                                break;
                        case 'G':
                        case 'g':
                                iseq->seq[i][j] = 2;
                                break;
                        case 'T':
                        case 't':
                                iseq->seq[i][j] = 3;
                                break;
                        default:
                                ERROR_MSG("Non ACGT letter in sequence:%d %s.",i,sequences[i]);
                                break;
                        }
                }
        }

	
        return iseq;
ERROR:
        free_ihmm_seq(iseq);
        return NULL;
}

void free_ihmm_seq(struct ihmm_sequences* iseq)
{
        if(iseq){
                free_2d((void**) iseq->seq);
                free_2d((void**) iseq->u);
                free_2d((void**) iseq->labels);
                MFREE(iseq->score);
                MFREE(iseq->len);
                MFREE(iseq);
        }
}


struct dp* init_dp(int num_states, int len)
{
        struct dp*  dp= NULL;
	
        MMALLOC(dp,sizeof(struct dp));
        dp->alloced_states = num_states;
        dp->alloced_len = len;
        dp->dp_matrix = NULL;
        dp->ID = NULL;
        dp->state_used = NULL;
	
	
        MMALLOC(dp->state_used, sizeof(int) * (dp->alloced_states));
        MMALLOC(dp->ID,sizeof(int) * (dp->alloced_states));
	
	
        dp->dp_matrix = malloc_2d_float(dp->dp_matrix , dp->alloced_len, dp->alloced_states, 0.0f);
	
	
        return dp;
ERROR:
        free_dp(dp);
        return NULL;
}

int clear_dp(struct dp* dp)
{
        int i,j;
	
        for(i = 0; i < dp->alloced_len;i++){
                for(j = 0; j < dp->alloced_states;j++ ){
                        dp->dp_matrix[i][j] = 0.0f;
                }
        }
        return OK;
}

int resize_dp(struct dp* dp, int num_states,int len)
{
        if(dp->alloced_states  < num_states || dp->alloced_len < len){
                dp->alloced_states = (num_states > dp->alloced_states)?  num_states+5 : dp->alloced_states  ;
                dp->alloced_len = (len > dp->alloced_len)?  len+10 : dp->alloced_len  ;
		
                MREALLOC(dp->state_used, sizeof(int) * (dp->alloced_states));
                MREALLOC(dp->ID,sizeof(int) * (dp->alloced_states));
		
                dp->dp_matrix = malloc_2d_float(dp->dp_matrix , dp->alloced_len, dp->alloced_states, 0.0f);
                return OK;
        }else{
                return OK;
        }
ERROR:
        free_dp(dp);
        return FAIL;
}

int free_dp(struct dp* dp)
{
        if(dp){
                free_2d((void**)dp->dp_matrix);
                MFREE(dp->state_used);
                MFREE(dp->ID);
                MFREE(dp);
        }
        return OK;
}


int select_random(float* vector,int num_elem,int* selected,float r )
{
        int i;
        float sum = 0.0f;
        //float r = rk_double(&model->rndstate);
        /*sum =0.0f;
          for(i = 0; i < num_elem;i++){
          sum += vector[i];
		
          //if(i +2 > num_elem){
		
          //DPRINTF3("%d %f: %f   %d",i, vector[i],sum, num_elem);
          //}
          }
          sum = fabsf((1.0f - sum)/sum);
          //FLT_EPSILON
	
          ASSERT(sum <= 0.00001 ,"Vector does not sum to 1:%20.20f",sum);
        */
        *selected = num_elem-1;
        sum =0.0f;
        for(i = 0; i < num_elem;i++){
		
                sum += vector[i];
                if(r <= sum){
                        *selected = i;
                        break;
                }
        }
        return OK;
//ERROR:
//	return FAIL;
}




#ifdef ITEST

static void pretty_print_vector_float(char* description, float* vector, int num,int log);
static void pretty_print_vector_int(char* description, int* vector, int num);
static int simulate_sequence_occasionally_dishonest_casino(char* seq, int len);

int main(const int argc,const char * argv[])
{
        fprintf(stdout,"Hello world\n");
	
        int i;
        //rk_randomseed( &rndstate);
        init_logsum();
	
        const char *tmp_seq[119] = {
                //CLUSTER:24832 (119)
                //	>Seq0_1
                "ACAGGCTAAAGGAGGGGGCAGTCCCCA",
                "AGGCTAAAGGAGGGGGCAGTCCCCACC",
                "AGGCTAAAGGAGGGGGCAGTCCCCACC",
                "AGTCCCCACCATATTTGAGTCTTTCTC",
                "AGTGGATATCACAGGCTAAAGGAGGGG",
                "AGTGGATATCACAGGCTAAAGGAGGGG",
                "AGTGGATATCACAGGCTAAAGGAGGGG",
                "AGTGGATATCACAGGCTAAAGGAGGGG",
                "AGTGGATATCACAGGCTAAAGGAGGGT",
                "CTAAAGGAGGGGGCAGTCCCCACCATA",
                "GAGGCTAAAGGAGGGGGCAGTCCCCAT",
                "GAGGGGGCAGTCCCCACCATATTTGAG",
                "GAGGGGGCAGTCCCCACCATATTTGAG",
                "GAGGGGGCAGTCCCCACCATATTTGAG",
                "GAGGGGGCAGTCCCCACCATATTTGAT",
                "GAGGGGGCAGTCCCCACCATATTTGAT",
                "GAGGGGGCAGTCCCCACCATATTTGAT",
                "GAGTCTTTCTCCAAGTTGCGCCGGACA",
                "GCAGTCCCCACCATATTTGAGTCTTTC",
                "GCAGTCCCCACCATATTTGAGTCTTTC",
                "GCAGTCCCCACCATATTTGAGTCTTTC",
                "GCAGTCCCCACCATATTTGAGTCTTTC",
                "GCAGTCCCCACCATATTTGAGTCTTTC",
                "GCAGTCCCCACCATATTTGAGTCTTTC",
                "GCAGTCCCCACCATATTTGAGTCTTTT",
                "GCAGTCCCCACCATATTTGAGTCTTTT",
                "GCAGTCCCCACCATATTTGAGTCTTTT",
                "GCAGTCCCCACCATATTTGAGTCTTTT",
                "GCAGTCCCCACCATATTTGAGTCTTTT",
                "GCAGTCCCCACCATATTTGAGTCTTTT",
                "GCAGTCCCCACCATATTTGAGTCTTTT",
                "GCAGTCCCCACCATATTTGAGTCTTTT",
                "GCAGTCCCCACCATATTTGAGTCTTTT",
                "GCAGTGGATATCACAGGCTAAAGGAGT",
                "GCTAAAGGAGGGGGCAGTCCCCACCAT",
                "GGAGGGGGCAGTCCCCACCATATTTGA",
                "GGAGGGGGCAGTCCCCACCATATTTGA",
                "GGATATCACAGGCTAAAGGAGGGGGCA",
                "GGCAGTCCCCACCATATTTGAGTCTTC",
                "GGCAGTCCCCACCATATTTGAGTCTTT",
                "GGCAGTCCCCACCATATTTGAGTCTTT",
                "GGCAGTCCCCACCATATTTGAGTCTTT",
                "GGCAGTCCCCACCATATTTGAGTCTTT",
                "GGCAGTCCCCACCATATTTGAGTCTTT",
                "GGCAGTCCCCACCATATTTGAGTCTTT",
                "GGCAGTCCCCACCATATTTGAGTCTTT",
                "GGCAGTCCCCACCATATTTGAGTCTTT",
                "GGCAGTCCCCACCATATTTGAGTCTTT",
                "GGCAGTCCCCACCATATTTGAGTCTTT",
                "GGCAGTCCCCACCATATTTGAGTCTTT",
                "GGCAGTCCCCACCATATTTGAGTCTTT",
                "GGCAGTCCCCACCATATTTGAGTCTTT",
                "GGCAGTCCCCACCATATTTGAGTCTTT",
                "GGGCAGTCACCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGGCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGGCAGTCCCCACCATATTTGAGTCT",
                "GGGGCAGTCCCCACCATATTTGAGTCT",
                "GGGGCAGTCCCCACCATATTTGAGTCT",
                "GGGGGCAGTCCCCACCATATTTGAGTC",
                "GGGGGCAGTCCCCACCATATTTGAGTC",
                "GTCCCCACCATATTTGAGTCTTTCTCT",
                "TCCCCACCATATTTGAGTCTTTCTCCA",
                "TCCCCACCATATTTGAGTCTTTCTCCA",
                "TGAGTCTTTCTCCAAGTTGCGCCGGAT",
                "TGGATATCACAGGCTAAAGGAGGGGGC"};
	
        char** sequences = NULL;
        //rk_randomseed( &rndstate);
        int len = 100;
        int numseq = 119;
	
        MCALLOC(sequences, numseq, sizeof(char*) );
        for(i = 0; i < numseq;i++){
                MMALLOC(sequences[i], sizeof(char)* (len+1));
	 
                snprintf(sequences[i], len, "%s",tmp_seq[i] );
        }
	
        struct iHMM_model* model = NULL;
	
        RUNP(model = init_iHMM_model());
        //hmm_data->max_len
        model->expected_K = 10;//  	  ; //iseq->max_len*2;//aln_len;//100;// iseq->len[0];
        model->alpha0 = 0.0;
        model->gamma = 0.0;
        model->alpha0_a  = 4.0;//(float) total_len;
        model->alpha0_b = 1.0;
	
	
        model->gamma_a = 3.0;// 3.0;//(float)  total_len * 0.05;
        model->gamma_b = 6.0;
	
	
        model->numb = 500;
        model->nums = 10;
        model->numi = 10;

        //DPRINTF2("START PGAS.");
        RUN(particle_gibbs_with_ancestors_controller(model,  sequences, numseq));
	
	
        print_emisson_matrix(model);
        
        //DPRINTF2("DONE PGAS.");
        free_iHMM_model(model);
	
        for(i = 0; i < numseq;i++){
                MFREE(sequences[i]);
        }

        MFREE(sequences);
	
        return EXIT_SUCCESS;
ERROR:
        return EXIT_FAILURE;
}
/*
  static int simulate_sequence_occasionally_dishonest_casino(char* seq, int len)
  {
	double transition[2][2] = { { 0.95, 0.05 }, { 0.1, 0.9}};
	
	double emission[2][4] =  {
  { 1.0/4.0, 1.0/4.0, 1.0/4.0, 1.0/4.0},
  { 0.1, 0.1, 0.1, 0.9}
	};
	//exit(0);
	int state = 0;
	int i,j;
	double r = 0;
	double sum = 0;
	
	
	
	int k = 2;
	int l = 4;
	for(i = 0; i < len;i++){
  r = rk_double(&rndstate);
  sum = 0.0;
  for(j = 0; j < k;j++){
  sum += transition[state][j];
  if(r < sum){
  state = j;
  //seq->true_labels[i] = state;
  break;
  }
  }
		
  r = rk_double(&rndstate);
  ///fprintf(stderr,"r:%f\n",r);
  sum = 0.0;
  for(j = 0; j < l;j++){
  sum += emission[state][j];
  //	fprintf(stderr,"r:%f	sum:%f\n",r,sum);
  if(r < sum){
  seq[i] = "ACGT"[j];
  break;
  }
  }
  //fprintf(stderr,"\n");
	}
	seq[len] = 0;
	return OK;
  }

*/
void pretty_print_vector_float(char* description, float* vector, int num,int log)
{
        int i;
        float sum = 0.0;
        if(!log){
                sum = 0.0;
                fprintf(stdout,"------------------------------\n%s\n",description);
                for(i = 0; i < num;i++){
                        sum +=vector[i];
                        fprintf(stdout," %0.3f", vector[i]);
                }
                fprintf(stdout,"		sum:%f\n",sum);
        }else{
                sum = prob2scaledprob(0.0);
                fprintf(stdout,"------------------------------\n%s\n",description);
                for(i = 0; i < num;i++){
                        sum = logsum(sum, vector[i]);
                        //sum +=vector[i];
                        fprintf(stdout," %0.3f", scaledprob2prob(vector[i]));
                }
                fprintf(stdout,"		sum:%f\n",scaledprob2prob(sum));
        }
        fflush(stdout);
}

void pretty_print_vector_int(char* description, int* vector, int num)
{
        int i;
        float sum = 0.0;
        fprintf(stdout,"------------------------------\n%s\n",description);
        for(i = 0; i < num;i++){
                sum +=vector[i];
                fprintf(stdout," %3d", vector[i]);
        }
        fprintf(stdout,"		sum:%f\n",sum);
        fflush(stdout);
}
#endif

