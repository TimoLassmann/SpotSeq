
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "tldevel.h"



int sanity_check_hmm(struct hmm*  hmm)
{
        int i,j;
	
        init_logsum();
        float sum = prob2scaledprob(0.0);
	
        for(j = 1; j < hmm->tindex[STARTSTATE][0];j++){
                if(hmm->train_transitions[STARTSTATE][hmm->tindex[STARTSTATE][j]]){
                        sum = logsum(sum, hmm->transitions[STARTSTATE][hmm->tindex[STARTSTATE][j]]);
                        //	fprintf(stdout,"%f ",scaledprob2prob(hmm->transitions[STARTSTATE][hmm->tindex[STARTSTATE][j]]));
                }
        }
	
        sum = (int) (10000.0 * scaledprob2prob(sum));
	
        //ASSERT(sum == 10000.0 ,"Start state transitions do not sum to 1:%20.20f",sum);
        for(i = 2; i < hmm->num_states;i++){
                sum = prob2scaledprob(0.0);
                for(j = 1; j < hmm->tindex[i][0];j++){
                        if(hmm->train_transitions[i][hmm->tindex[i][j]]){
                                //DPRINTF2("%f	%d -> %d\n",hmm->transitions[i][hmm->tindex[i][j]],i,hmm->tindex[i][j] );
                                sum = logsum(sum, hmm->transitions[i][hmm->tindex[i][j]]);
                        }
                }
                sum = (int) (10000.0 * scaledprob2prob(sum));
                //ASSERT(sum == 10000.0 ,"State transitions from %d do not sum to 1:%20.20f",i,sum);
        }
	
        for(i = 2; i < hmm->num_states;i++){
                sum = prob2scaledprob(0.0);
                for(j = 0; j < 4;j++){
                        sum = logsum(sum, hmm->emissions[i][j]);
                }
		
                sum = (int) (10000.0 * scaledprob2prob(sum));
                //ASSERT(sum == 10000.0 ,"State transitions from %d do not sum to 1:%20.20f",i,sum);
        }
        return OK;
}


struct hmm* malloc_hmm(int num_states, int alphabet_len, int max_seq_len)
{
        struct hmm* hmm = NULL;
        int i = 0;
        int j = 0;
        MMALLOC(hmm, sizeof(struct hmm));
        hmm->alphabet_len =  alphabet_len;
        hmm->num_states = num_states + NUM_ADDITIONAL_STATES;
        hmm->max_seq_len = max_seq_len;
        hmm->emissions = NULL;
        hmm->transitions = NULL;
        hmm->tindex = NULL;
	
        hmm->emissions_e = NULL;
        hmm->transitions_e = NULL;
	
        //hmm->data = NULL;
        hmm->background = NULL;
	
        hmm->F = NULL;
        hmm->B = NULL;
	
        //hmm->F_memory = NULL;
        //hmm->B_memory = NULL;
	
        hmm->train_emissions = NULL;
        hmm->train_transitions = NULL;
	
        MMALLOC(hmm->train_emissions , sizeof(char) * hmm->num_states);
	
        MMALLOC(hmm->train_transitions,sizeof(char*) * hmm->num_states);
	
        //MMALLOC(hmm->emissions,sizeof(float*) * num_states);
        //MMALLOC(hmm->transitions,sizeof(float*) * num_states);
        MMALLOC(hmm->tindex,sizeof(float*) * hmm->num_states);
        MMALLOC(hmm->emissions_e,sizeof(float*) * hmm->num_states);
        MMALLOC(hmm->transitions_e,sizeof(float*) * hmm->num_states);
	
        MCALLOC(hmm->background, hmm->alphabet_len,sizeof(float));
	
	
        hmm->emissions = malloc_2d_float(hmm->emissions, hmm->num_states, hmm->alphabet_len, prob2scaledprob( 0.0f));
        hmm->transitions = malloc_2d_float(hmm->transitions, hmm->num_states, hmm->num_states,  prob2scaledprob(0.0f));
	
	
        for(i = 0; i < hmm->num_states;i++){
                hmm->train_emissions[i] = 1;
                //hmm->emissions[i] = NULL;
                //hmm->transitions[i] = NULL;
                hmm->emissions_e[i] = NULL;
                hmm->transitions_e[i] = NULL;
                hmm->tindex[i] = NULL;
                hmm->train_transitions[i] = NULL;
                MMALLOC(hmm->train_transitions[i],sizeof(char) * hmm->num_states);
                //hmm->F[i] = NULL;
                //hmm->B[i] = NULL;
		
                //MMALLOC(hmm->transitions[i],sizeof(float) * num_states);
                MMALLOC(hmm->transitions_e[i],sizeof(float) * hmm->num_states);
                MMALLOC(hmm->tindex[i],sizeof(float*) * (hmm->num_states+1));
		
		
                //MMALLOC(hmm->emissions[i], sizeof(float) * hmm->alphabet_len);
                MMALLOC(hmm->emissions_e[i], sizeof(float) * hmm->alphabet_len);
		
                for(j = 0; j < hmm->num_states;j++){
                        hmm->train_transitions[i][j] = 1;
                }
                //MMALLOC(hmm->F[i],sizeof(float) * max_seq_len);
                //MMALLOC(hmm->B[i],sizeof(float) * max_seq_len);
		
		
        }
        hmm->F = malloc_2d_float(hmm->F,max_seq_len, hmm->num_states,  0.0);
        hmm->B = malloc_2d_float(hmm->B,max_seq_len, hmm->num_states , 0.0);
	
        return hmm;
ERROR:
        if(hmm){
                for(i = 0; i < hmm->num_states;i++){
                        MFREE(hmm->train_transitions[i]);// ,sizeof(char) * num_states);
                        MFREE(hmm->transitions_e[i]);// ,sizeof(float) * num_states);
                        MFREE(hmm->tindex[i]);// ,sizeof(float*) * (num_states+1));
                        MFREE(hmm->emissions_e[i]);// , sizeof(float) * hmm->alphabet_len);
                }
                if(hmm->emissions){
                        free_2d((void**) hmm->emissions );
                }
                if(hmm->transitions){
                        free_2d((void**) hmm->transitions );
                }
		
                free_2d((void**) hmm->F);
                free_2d((void**) hmm->B);
		
		
                MFREE(hmm->train_emissions);// , , sizeof(char) * num_states);
		
                MFREE(hmm->train_transitions);// ,,sizeof(char*) * num_states);
		
                //MMALLOC(hmm->emissions,sizeof(float*) * num_states);
                //MMALLOC(hmm->transitions,sizeof(float*) * num_states);
                MFREE(hmm->tindex);// ,,sizeof(float*) * num_states);
                MFREE(hmm->emissions_e);// ,,sizeof(float*) * num_states);
                MFREE(hmm->transitions_e);// ,,sizeof(float*) * num_states);
		
                MFREE(hmm->background);// ,, sizeof(float) *hmm->alphabet_len);
		
                MFREE(hmm);// ,, sizeof(struct hmm));
        }
        return NULL;
}

void free_hmm(struct hmm* hmm)
{
        int i;
	
	
	
	
        for(i = 0; i < hmm->num_states;i++){
                if(hmm->emissions[i]){
                        //MFREE(hmm->emissions[i]);
                        MFREE(hmm->emissions_e[i]);
                }
                MFREE(hmm->train_transitions[i]);
                //MFREE(hmm->transitions[i]);
                MFREE(hmm->transitions_e[i]);
                MFREE(hmm->tindex[i]);
		
        }
        free_2d((void**)hmm->emissions );
        free_2d((void**)hmm->transitions );
        MFREE(hmm->train_transitions);
        MFREE(hmm->train_emissions);
	
        free_2d((void**) hmm->F);
        free_2d((void**) hmm->B);
	
        //MFREE(hmm->F);
        //MFREE(hmm->B);
        //MFREE(hmm->F_memory);
        //MFREE(hmm->B_memory);
	
        //MFREE(hmm->transitions);
        MFREE(hmm->tindex);
        MFREE(hmm->transitions_e);
        //MFREE(hmm->emissions);
        MFREE(hmm->emissions_e);
        MFREE(hmm->background);
        MFREE(hmm);
}

struct hmm* forward(struct hmm* hmm, char* a, int len)
{
        int i,j,c,f;
	
        float** matrix = hmm->F;
        float* last= 0;
        float* cur = 0;
        const float* trans = 0;
	
        float tmp = 0;
	
        cur = matrix[0];
	
        //for(i = 0; i < hmm->num_states;i++){
        for(j = 0; j < hmm->num_states;j++){
                cur[j]  = -INFINITY;
        }
        cur[STARTSTATE] = 0.0f;
	
        for(i = 1; i < len+1;i++){
                last = cur;
                cur = matrix[i];
                for(j = 0; j < hmm->num_states;j++){
                        cur[j] = -INFINITY;
                }
                for(j = 0; j < hmm->num_states;j++){
                        tmp = last[j];
                        trans = hmm->transitions[j];
                        for(c = 1; c < hmm->tindex[j][0];c++){
                                f = hmm->tindex[j][c];
                                cur[f] = logsum(cur[f], tmp + trans[f] );//+ hmm->emissions[c][(int)a[i-1]]);
                        }
                }
                for(c = 2;c < hmm->num_states;c++){
                        cur[c] += hmm->emissions[c][(int)a[i-1]];
                }
        }
        //All goes to 1.
        last = cur;//matrix[len];
        cur = matrix[len+1];
	
	
        for(j = 0; j < hmm->num_states;j++){
                cur[j] = -INFINITY;// prob2scaledprob(0.0);
        }
	
	
        for(j = 2; j < hmm->num_states;j++){
                cur[ENDSTATE] = logsum(cur[ENDSTATE],last[j] + hmm->transitions[j][ENDSTATE]);
        }
        hmm->f_score = cur[ENDSTATE];// matrix[ENDSTATE][i];
        return hmm;
}


struct hmm* backward(struct hmm* hmm, char* a, int len)
{
        int i,j,c,f;
	
        float** matrix = hmm->B;
	
        float* next= 0;
        float* cur = 0;
        const float* trans = 0;
	
        cur = matrix[len+1];
	
        for(j = 0; j < hmm->num_states;j++){
                cur[j] = -INFINITY;
        }
	
        cur[ENDSTATE] = 0.0f;
	
        next = cur;
	
        cur = matrix[len];
        for(j = 0; j < hmm->num_states;j++){
                cur[j] =  hmm->transitions[j][ENDSTATE] + next[ENDSTATE];
        }
        for(c = 2;c < hmm->num_states;c++){
                cur[c] += hmm->emissions[c][(int)a[len-1]];
        }
	
        // backward recursion...
        for(i = len-1; i > 0; i -- ){
                next = cur;
                cur = matrix[i];
                for(j = 0; j < hmm->num_states;j++){
                        trans = hmm->transitions[j];
                        cur[j] = -INFINITY;
                        for(c = 1; c < hmm->tindex[j][0];c++){
                                f = hmm->tindex[j][c];
                                cur[j] = logsum(cur[j],trans[f] + next[f]);// hmm->emissions[c][(int)a[i]]);// + next[c]);
                        }
                }
                for(j = 2; j < hmm->num_states;j++){
                        cur[j] += hmm->emissions[j][(int)a[i-1]];
                }
        }
	
        cur = matrix[0];
        next = matrix[1];
	
        for(j = 0; j < hmm->num_states;j++){
                cur[j] = -INFINITY;// prob2scaledprob(0.0f);
        }
        for(i = 0; i < hmm->num_states;i++){
                cur[STARTSTATE] = logsum(cur[STARTSTATE], hmm->transitions[STARTSTATE][i] + next[i]);//  + hmm->emissions[i][(int)a[0]]  );
        }
        hmm->b_score = cur[STARTSTATE];
        return hmm;
}


struct hmm* posterior_decoding(struct hmm* hmm, char* a, int len,int* path)
{
        int i,j,c,f,best;
        float** Fmatrix = hmm->F;
        float** Bmatrix = hmm->B;
	
        float* this_F = 0;
        float* this_B = 0;
	
        //const float* trans = 0;
        float max = prob2scaledprob(0.0);
        float total = hmm->f_score;
	
        for(i = 1; i <= len+1;i++){
                //last_F = Fmatrix[i-1];
                this_F = Fmatrix[i];
                this_B = Bmatrix[i];
                for(j = 0; j < hmm->num_states;j++){
                        this_F[j] = scaledprob2prob((this_F[j]  +( this_B[j] -hmm->emissions[j][(int)a[i-1]] ))   -total);
                }
        }
	
        for(j = 0; j < hmm->num_states;j++){
                Fmatrix[len][j] = Fmatrix[len][j]  + Fmatrix[len+1][ENDSTATE] ;
                Bmatrix[len][j] = ENDSTATE;
        }
        best = -1;
        /* Fill B with best transition pointers... */
        for(i = len-1; i >= 0; i-- ){
                for(j = 0; j < hmm->num_states;j++){
                        //trans = hmm->transitions[j];
                        max = -INFINITY;
                        for(c = 1; c < hmm->tindex[j][0];c++){
                                f = hmm->tindex[j][c];
                                if( Fmatrix[i+1][f] > max){
                                        max  = Fmatrix[i+1][f] ;
                                        best = f;
                                }
                        }
                        Fmatrix[i][j] = max + Fmatrix[i][j];
                        Bmatrix[i][j] = best;
                }
        }
       
        int state = STARTSTATE;
	
	
        /* traceback */
        i = 0;
        while (i <= len){
                //going to
                if(i){
                        path[i-1] = state;
                }
		
                state = hmm->B[i][state];
                i++;
        }
	

        for(i =0; i < len;i++){
                fprintf(stdout,"%2d %2d %2d\n", i,(int) a[i], path[i] );
        }
        fprintf(stdout,"\n");
	
	
	
	
        return hmm;
}
 
