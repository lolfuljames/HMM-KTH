import java.util.*;
import java.io.*;

public class HMM3{
    
    private static Scanner sc = new Scanner(System.in);
    public static void main(String[] args){

        final int MAX_ITERATIONS = 5000; // maximum of 500,000 iterations
        final double EPSILON = 1e-7;
       
        /*
         * Matrix where row = current state, col = next state (timestep), elements = joint probability of state entering next state (row -> col).
         */
        Matrix transition_matrix = create_matrix();

        /*
         *  Matrix where row = current state, col = observation, elements = joint probability of state and observation.
         */
        Matrix emission_matrix = create_matrix();

        /*
         *  Matrix where row size = 1, col = states, elements = probability of entering initial state (col).
         */
        Matrix initial_state_matrix = create_matrix();
        
        /*
         *  Array of observations where array[0] is the first observation, array[n-1] is the n-th observation
         */
        int[] observation_sequence = create_observation_sequence();

        /*
         *  Matrix where row = alpha at time(row), col = states, elements = probability of entering state (col) given sequence of previous observations.
         */
        ArrayList<Matrix> alpha_matrix = new ArrayList<Matrix>();

        /*
         * Array of double that represent the scale factor used. Where scale factor is scaling from the previous timestep.
         */
        double[] scale_factor = new double[observation_sequence.length];

        /*
         *  Matrix where row = beta at time(row), col = states, elements = probability of observing future observations at current timestep(row) and state(col).
         */
        ArrayList<Matrix> beta_matrix = new ArrayList<Matrix>();

        /*
         *  Array of matrix where index = timestep, where row = current state, col = future state, elements = probability of being in state row, and next state in timestep t+1 being state col.
         */
        ArrayList<Matrix> digamma_matrix = new ArrayList<Matrix>();

        /*
         *  Array of matrix where index = timestep, where col = current state, elements = probability of being in state row.
         */
        ArrayList<Matrix> gamma_matrix = new ArrayList<Matrix>();

        double numer = 0;
        double denom = 0;
        double logprob = -10000000;
        double old_logprob = -100000000;
        int iteration = 0;

        // System.out.println(transition_matrix);
        // System.out.println(emission_matrix);
        // System.out.println(observation_sequence);

        while(iteration<MAX_ITERATIONS && logprob>old_logprob+EPSILON){
//            System.out.println(iteration);
            old_logprob = logprob;


            alpha_pass(alpha_matrix, observation_sequence, transition_matrix, emission_matrix, initial_state_matrix, scale_factor);
            beta_pass(beta_matrix, observation_sequence, transition_matrix, emission_matrix, initial_state_matrix, scale_factor);
            // System.out.println(beta_matrix.get(beta_matrix.size()-1));
            // System.out.println(alpha_matrix.get(alpha_matrix.size()-1));
            // System.out.println(scale_factor[0]);
            // sc.nextLine();
            di_gamma(digamma_matrix, gamma_matrix, alpha_matrix, beta_matrix, observation_sequence, transition_matrix, emission_matrix);
            // System.out.println(digamma_matrix.get(0));
            // System.out.println(digamma_matrix.get(digamma_matrix.size()-1));

            Matrix final_alpha = alpha_matrix.get(alpha_matrix.size()-1);
            double total = 0;

            // for(int i=0; i<final_alpha.get_elements()[0].length; i++){
            //     total += final_alpha.get_elements()[0][i];
            // }
            // for(int i=0; i<scale_factor.length; i++){
            //     total *= scale_factor[i];
            // }

            // System.out.println(total);

            //re-estimating initial state matrix
            initial_state_matrix = gamma_matrix.get(0);

            //re-estimating transition matrix
            for(int i=0; i< transition_matrix.get_row(); i++){
                for(int j=0; j< transition_matrix.get_col(); j++){
                    numer = 0;
                    denom = 0;
                    for (int time=0; time< observation_sequence.length-1; time++){
                        numer += digamma_matrix.get(time).get_elements()[i][j];
                        denom += gamma_matrix.get(time).get_elements()[0][i];
                    }
                    transition_matrix.set_element(numer/denom, i, j);
                }
            }

            //re-estimating emission matrix
            for(int i=0; i< emission_matrix.get_row(); i++){
                for(int j=0; j< emission_matrix.get_col(); j++){
                    numer = 0;
                    denom = 0;
                    for (int time=0; time< observation_sequence.length; time++){
                        if(j==observation_sequence[time]) numer += gamma_matrix.get(time).get_elements()[0][i];
                        denom += gamma_matrix.get(time).get_elements()[0][i];
                    }
                    emission_matrix.set_element(numer/denom, i, j);
                }
            }

            //calculating probability
            logprob = 0;
            for (int i =0; i<observation_sequence.length; i++){
                logprob += Math.log(1/scale_factor[i]);
            }

//            System.out.println(digamma_matrix.get(0));
            logprob *= -1;
            iteration ++;
            alpha_matrix = new ArrayList<Matrix>();
            beta_matrix = new ArrayList<Matrix>();
            digamma_matrix = new ArrayList<Matrix>();
            gamma_matrix = new ArrayList<Matrix>();
            scale_factor = new double[observation_sequence.length];
            // sc.nextLine();
        }
//        System.out.println(digamma_matrix.get(0));
         System.out.println(transition_matrix);
         System.out.println(emission_matrix);
        // System.out.println(logprob);
        // System.out.println(old_logprob);
    }

    public static Matrix create_matrix(){
        int row;
        int col;

        row = sc.nextInt();
        col = sc.nextInt();
        double[][] elements = new double[row][col];

        for(int j=0; j<row; j++){
            for(int k=0; k<col; k++){
                elements[j][k] = sc.nextDouble();
            }
        }
        Matrix new_matrix = new Matrix(row,col,elements);
        return new_matrix;
    }

    public static int[] create_observation_sequence(){
        int[] observation_sequence;
        int length = sc.nextInt();

        observation_sequence = new int[length];
        for(int i=0; i<length; i++){
            int result = sc.nextInt();
            observation_sequence[i] = result;
        }

        return observation_sequence;
    }

    /*
     * Function used to initialise forward algorithm
     */
    public static void alpha_pass(ArrayList<Matrix> alpha_matrix, int[] observation_sequence, Matrix transition_matrix, Matrix emission_matrix, Matrix initial_state_matrix, double[] scale_factor){
        
        Matrix next_state = new Matrix();
        Matrix current_state = new Matrix();
        double current_scale = 0;
        
        while(alpha_matrix.size() != observation_sequence.length){
            //calculate alpha(1)
            // alpha(1) = pi (written below) *P(x|0)
            if(alpha_matrix.size() == 0) next_state = initial_state_matrix;//.multiply(transition_matrix);
            // alpha(t) = alpha(t-1)*B (written below) *P(x|0)
            else next_state = current_state.multiply(transition_matrix);

            // multiplication of P(x|0)
            next_state = alpha_multiply_observation(next_state, emission_matrix, observation_sequence[alpha_matrix.size()]);

            // Obtaining summation of alpha(t) to get scale factor
            current_scale = 0;
            for(int i = 0; i<next_state.get_col(); i++){
                current_scale += next_state.get_elements()[0][i];
            }

            // Scaling by scale factor
            for(int i = 0; i<next_state.get_col(); i++){
                next_state.get_elements()[0][i]/=current_scale;
            }

            // Storing alpha and scale factor
            scale_factor[alpha_matrix.size()] = current_scale;
            alpha_matrix.add(next_state);

            // alpha is passed down and saved
            current_state = next_state;
        }
    }

    /*
     * Function used to multiply given state distribution probablility and emission distribution of states, similar to dot product of emissions.
     */

    public static Matrix alpha_multiply_observation(Matrix state_probability, Matrix emission_matrix, int observation){

        double[][] alpha = new double[1][emission_matrix.get_row()];
        for(int i=0; i<emission_matrix.get_row(); i++){
            for(int j=0; j<emission_matrix.get_col(); j++){
                if(observation == j){
                    alpha[0][i] = Math.exp(Math.log(state_probability.get_elements()[0][i]) + Math.log(emission_matrix.get_elements()[i][j]));
                } 
            }
        }
        return  new Matrix(alpha.length,alpha[0].length,alpha);
    }

    public static void beta_pass(ArrayList<Matrix> beta_matrix, int[] observation_sequence, Matrix transition_matrix, Matrix emission_matrix, Matrix initial_state_matrix, double[] scale_factor){

        //initializing beta(T), current is t, next is t+1
        Matrix current_beta;
        Matrix next_beta;

        double[][] final_beta = new double[1][emission_matrix.get_row()];
        for(int i = 0; i<final_beta[0].length; i++){
            final_beta[0][i] = 1/scale_factor[observation_sequence.length-1];
        }
        current_beta = new Matrix(final_beta.length, final_beta[0].length, final_beta);

        //initializing arraylist of beta matrices
        for(int i =0; i<observation_sequence.length; i++){
            beta_matrix.add(new Matrix(1,transition_matrix.get_row()));
        }

        beta_matrix.set(observation_sequence.length-1, current_beta);
        // System.out.println("Timestep 0 Beta:");
        // System.out.println(current_beta);

        // for(int i = observation_sequence.length-2; i>= 0; i--){
        //     next_beta = current_beta;
        //     current_beta = next_beta.multiply(transition_matrix);
        //     current_beta = alpha_multiply_observation(current_beta, emission_matrix, observation_sequence[i+1]);
        //     current_beta = current_beta.multiply(1/scale_factor[i]);
        //     beta_matrix.set(i, current_beta);
        // }
        
        double result = 0;
        for(int time = observation_sequence.length-2; time>= 0; time--){
            for(int i = 0; i<transition_matrix.get_row(); i++){
                current_beta = beta_matrix.get(time);
                result = 0;
                for(int j=0; j< transition_matrix.get_col(); j++){
                    result = result + transition_matrix.get_elements()[i][j]*
                        emission_matrix.get_elements()[j][observation_sequence[time+1]]*
                        beta_matrix.get(time+1).get_elements()[0][j];
                }
                result = result / scale_factor[time];
                current_beta.set_element(result, 0, i);
                beta_matrix.set(time, current_beta);
            }
        }

    }

    public static void di_gamma(ArrayList<Matrix> digamma_matrix, ArrayList<Matrix> gamma_matrix, ArrayList<Matrix> alpha_matrix,ArrayList<Matrix> beta_matrix,int[] observation_sequence,Matrix transition_matrix,Matrix emission_matrix){
        double[][] di_gamma = new double[transition_matrix.get_row()][transition_matrix.get_col()];
        double[][] gamma = new double[1][transition_matrix.get_row()];
        double denom;

//        System.out.println("ALPHA:::");
//        System.out.println(alpha_matrix.get(0));
//        System.out.println("BETA:::");
//        System.out.println("Timestep 0 Beta:");
//        System.out.println(beta_matrix.get(0));

        for(int time=0; time<observation_sequence.length-1; time++){
            denom =0;
            double a,b,c,d,e;
            di_gamma = new double[transition_matrix.get_row()][transition_matrix.get_row()];
            gamma = new double[1][transition_matrix.get_row()];
            for(int i=0; i<transition_matrix.get_row(); i++){
                for(int j=0; j<transition_matrix.get_col(); j++){
//                    a=alpha_matrix.get(time).get_elements()[0][i];
//                    b=transition_matrix.get_elements()[i][j];
//                    c=emission_matrix.get_elements()[j][observation_sequence[time+1]];
//                    d=beta_matrix.get(time+1).get_elements()[0][j];
//                    e=a*b*c*d;
                    denom = denom + alpha_matrix.get(time).get_elements()[0][i]*transition_matrix.get_elements()[i][j]*
                        emission_matrix.get_elements()[j][observation_sequence[time+1]]*beta_matrix.get(time+1).get_elements()[0][j];
                }
            }
//            if(time ==998) System.out.printf("DENOM HERE %f", denom);
            for(int i =0; i<transition_matrix.get_row(); i++){
                for(int j=0; j<transition_matrix.get_col(); j++){
//                    a=alpha_matrix.get(time).get_elements()[0][i];
//                    b=transition_matrix.get_elements()[i][j];
//                    c=emission_matrix.get_elements()[j][observation_sequence[time+1]];
//                    d=beta_matrix.get(time+1).get_elements()[0][j];
//                    e=a*b*c*d;
                    di_gamma[i][j] = alpha_matrix.get(time).get_elements()[0][i]*transition_matrix.get_elements()[i][j]*
                        emission_matrix.get_elements()[j][observation_sequence[time+1]]*
                        beta_matrix.get(time+1).get_elements()[0][j]/denom;
                    gamma[0][i] += di_gamma[i][j];
                }
            }
            digamma_matrix.add(new Matrix(di_gamma.length,di_gamma[0].length,di_gamma));
            gamma_matrix.add(new Matrix(gamma.length,gamma[0].length,gamma));
        }

        //special case for gamma(t-1)
        denom = 0;
        for(int i =0; i<emission_matrix.get_row(); i++){
            denom += alpha_matrix.get(alpha_matrix.size()-1).get_elements()[0][i];
        }
        gamma = new double[1][emission_matrix.get_row()];
        for(int i =0; i<emission_matrix.get_row(); i++){
            gamma[0][i] = alpha_matrix.get(alpha_matrix.size()-1).get_elements()[0][i]/denom;
        }
        gamma_matrix.add(new Matrix(gamma.length,gamma[0].length,gamma));
    }
}