/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package sunflecks;
import au.com.bytecode.opencsv.CSVReader;
import au.com.bytecode.opencsv.CSVWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;
import static java.lang.Math.*;
import javax.swing.*;
import org.math.plot.*;

/**
 *
 * @author lfhan
 */
public final class Sunflecks {

    //stomata para
    public double gsmax = 0.05;
    public double gsmin = 0.007;
    public double gb = 1440000;
    public double alphag = 1.0928E-03;
    public double thetag = 9.9745E-01;
    /**/
    public double rd = 2; //miumol
    public double vjmax = 24.2; //miumol  = 2 * this.vcmax;
    public double vcmax = 76.7764718311057; //miumol
    public double vfmax = 96.5000200233892; //this.vfmax = 1.8 * this.vcmax;
    public double vfmin = 0.05;
    public double vcmin = 0.01;
    /**/
    public double alphac = 0.00339868967379728; //est.
    public double thetac = 0.884412945823268;
    public double alphaf = 0.3;
    public double thetaf = 0.95;
    public double alphaj = 0.00542924797508813;
    public double thetaj = 0.074110803049945;
    /**/
    public double rmax = 370; //miu mol, also used in dy function
    public double tmax = 400; //miu mol, also used in dy function
    public double kr = 5; //need est. miu mol, also used in dy function
    public double kt = 5;  //need est. miu mol, also used in dy function 
    public double gamma = 4.44;//also used in dy function
    public double ka = 1;
    public double kc = 31;
    public double ko = 15500;
    public double patm = 100000;
    public double ca = 39;//changed
    public double po2 = 21000;
    /**/
    public double taugi = 215.5893;
    public double taugd = 324.2205;
    public double taupi = 206.6032;
    public double tauw = 208.7749;
    
    //psyn timeconst
    public double psi = 0.321430269302047;
    double tauci = 79.8627652635799;
    double taucd = 110.598886780825;
    double taufi=471.475996185132;
    double taufd=304.671929390908;
    
    
    //State variables
    public double signal;
    public double pi;
    public double water;
    
    public double poolT;
    public double poolR;
    public double poolG;
        
    public double ci;
    public double irr;
    public double gs;
    public double vj;
    public double vf;
    public double vc;
    public double ass;
    public double seq;

    
    


    public static void run(double iniPAR) throws FileNotFoundException, IOException {
        //instantiate the class, initialize it with equilibrium at PAR=100
        Sunflecks simSunflecks=new Sunflecks(iniPAR);       
        double dt=0.01;        //time step

        CSVReader reader = new CSVReader(new FileReader("input.csv"));
        CSVWriter writer = new CSVWriter(new FileWriter("output.csv"), '\t');
                
        //read time and PAR from data file and put into two arrays
        List<String[]> allElements = reader.readAll();
        int nDatapts=allElements.size();
        double[] timearr=new double[nDatapts];
        double[] irrarr=new double[nDatapts];
        double[] assarr=new double[nDatapts];
        String[][] tmparr = allElements.toArray(new String[0][]);    
        for (int i=0;i<nDatapts;i++){
            timearr[i]=Double.parseDouble(tmparr[i][0]);
            irrarr[i]=Double.parseDouble(tmparr[i][1]);        
        }
               
        
        
        double curtime=0;
        for (int i=0;i<nDatapts;i++){
            while (curtime<timearr[i]){
                //do the whole simulation every time step
                curtime+=dt; 
                simSunflecks.calcDyn(irrarr[i]);             
            }
            //get An at each data point and write into output file
            assarr[i]=simSunflecks.ass;
            String str2 =  String.valueOf(simSunflecks.ass);
            String str1 =  String.valueOf(timearr[i]);
            String[] strarr= {str1,str2};
            writer.writeNext(strarr);
        }        
        writer.close();
        double teststr=simSunflecks.ass;
        System.out.println(curtime);
        System.out.println(teststr);
        
        
        // create your PlotPanel (you can use it as a JPanel)
        Plot2DPanel plot = new Plot2DPanel();
        // define the legend position
        plot.addLegend("SOUTH");
        // add a line plot to the PlotPanel
        //plot.addLinePlot("my plot", timearr, irrarr);
        plot.addLinePlot("another plot", timearr, assarr);
        // put the PlotPanel in a JFrame like a JPanel
        JFrame frame = new JFrame("a plot panel");
        frame.setSize(600, 600);
        frame.setContentPane(plot);
        frame.setVisible(true);
    }
    
    //constructor be called when instantiate the class
    public Sunflecks(double iniIRR){
        //Assume the leaf is in equilibrium at the initial PAR 
        //initialize it by calling steady state model
        this.calcSteady(iniIRR);
    }
    
    public void calcSteady(double irr){
        this.irr=irr;
        this.vj=getvxeq(this.alphaj,this.thetaj,this.vjmax,0);
        this.vf=getvxeq(this.alphaf,this.thetaf,this.vfmax,this.vfmin);
        this.vc=getvxeq(this.alphac,this.thetac,this.vcmax,this.vcmin);
        updateStdygs();
        double newci=30; //first guess
        int ici=0;
        while(abs(newci-this.ci) > 0.0001 && ici<1000) {
            ici = ici + 1;
            if (ici==999){
                System.out.println("ici error");
            }
            //double oldci=this.ci;
            this.ci= newci;
            double fx1=fci(this.ci);
            double fx2=fci(this.ci-0.01);
            //Netwon's Method: xnew=xold-f(xold)/f'(xold)
            newci=this.ci-fx1*(0.01/(fx1-fx2));        
        }
        double gg = getgg(); 
        this.ass=(this.ca-this.ci)*gg/this.patm;                 
    }
    
    public void calcDyn(double irr){
        this.irr=irr;
        double dt=0.01;
        updateDyngs(dt);
        double gg=getgg();
        this.ci = this.ca - ((this.ass * this.patm) / gg);
        double newci=this.ci+1; //first guess
        for(int ci_i=1; ci_i<500; ci_i++){
            double oldci=this.ci;
            this.ci=newci;
            double fx1=dynfci(this.ci);
            double fx2=dynfci(oldci);
            //secent method
            newci=this.ci-fx1*((this.ci-oldci)/(fx1-fx2));
            if (abs(newci - this.ci) < 0.0001) {
                break;
            }
        }
        if (abs(newci - this.ci) > 0.0001) {
            System.out.println("dyn ci error");
        }
        this.vj = getvxeq(this.alphaj, this.thetaj, this.vjmax, 0);
        double vfeq=getvxeq(this.alphaf,this.thetaf,this.vfmax,this.vfmin);
        double vceq=getvxeq(this.alphac,this.thetac,this.vcmax,this.vcmin);
        double wc=getwc(this.ci);
        double[] tau={this.tauci,this.taucd,this.taufi,this.taufd}; //time const
        //create a object of Differential Equation
        Dfun photoDif=new Dfun(tau,this.rmax, this.tmax, this.kr,
                this.kt, this.gamma, this.psi,this.vj,vfeq,vceq, this.ci, wc);
        //initial value of state variables
        double[] y0={this.vc, this.vf, this.poolT, this.poolR, this.poolG};
        //calculate those values after one time step
        double[] yf=rk4(y0,dt,photoDif);
        this.vc=yf[0];
        this.vf=yf[1];
        this.poolT=yf[2];
        this.poolR=yf[3];
        this.poolG=yf[4];
        wc=getwc(this.ci);
        this.ass=getass1(wc); 
        System.out.println(this.ass);
    }
    
    //dynamic function of ci which we want it to be zero
    public double dynfci(double localci){
        double wc=getwc(localci);
        double ass1 = getass1(wc);
        double gg = getgg();
        double ass2=(this.ca-localci)*gg/this.patm;
        return ass1-ass2;
    }
    
    //steady state function of ci which we want it to be zero
    public double fci(double localci){
        double wc = getwc(localci);
        this.poolT=this.tmax;
        this.poolR = getR(); 
        double eqn2 = geteqn2(localci,wc);    //eqn2 in Pearcy1997 
        //bisection method 
        //first find T1 and T2 which statify eqn2(T1)<0 , eqn2(T2)>0 
        double T2 = 0;
        double T1 = 0;
        if (eqn2>0){
            T2=this.poolT;
            //find t1
            for(int findt1_i=1; findt1_i<500; findt1_i++){
                double Tportion=T2/findt1_i+1;             
                for(int findt1_j=1; findt1_j<findt1_i+1; findt1_j++){
                    T1=T2-Tportion*findt1_j;
                    this.poolT=T1;
                    this.poolR=getR();
                    eqn2 = geteqn2(localci,wc);
                    if (eqn2<0) {break;}                                        
                }
                if (eqn2<0) {break;}
            }
            if (eqn2>0){
                System.out.println("eqn2 error");
            }              
        }else{ //eqn2<0
            T1=this.poolT;
            //find t2
            for(int findt2_i=1; findt2_i<500; findt2_i++){
                double Tportion=T1/findt2_i+1;             
                for(int findt2_j=1; findt2_j<findt2_i+1; findt2_j++){
                    T2=T1-Tportion*findt2_j;
                    this.poolT=T2;
                    this.poolR=getR();
                    eqn2 = geteqn2(localci,wc); 
                    if (eqn2>0) {break;}                                        
                }
                if (eqn2>0) {break;}
            }
            if (eqn2<0){
                System.out.println("eqn2 error");
            }    
        }
        
        //bisection
        for (int bisec_i=1; bisec_i<1000; bisec_i++){
            if (abs(eqn2)<0.001){break;}
            this.poolT = (T1 + T2) / 2;
            this.poolR=getR();
            eqn2=geteqn2(localci,wc); 
            if (eqn2>0){
                T2=this.poolT;}
            else{
                T1=this.poolT;
            }   
        }
        if(abs(eqn2)>0.001){
             System.out.println("bisec error");
        }
        this.poolG = ((2 * this.gamma / localci) * wc 
                * (this.poolR / (this.kr + this.poolR))) / this.psi;  
        double ass1 = getass1(wc);
        double gg = getgg();
        double ass2=(this.ca-localci)*gg/this.patm;
        return ass1-ass2;
    }
    
    //method to update state varibles of steady state stomata
    public void updateStdygs(){
        this.seq=getSeq();
        this.gs=this.seq*this.gsmax;
        this.signal=this.seq;
        this.pi=this.seq;
        this.water=this.seq;
    }
    
    //method to update state varibles of dynamic stomata
    public void updateDyngs(double dt){
        this.seq=getSeq();
        Dfun gsDif=new Dfun(this.taugi,this.taugd,this.taupi,this.tauw,this.seq);
        double[] y0={this.signal,this.pi,this.water};
        double[] yf=rk4(y0,dt,gsDif);
        this.signal=yf[0];
        this.pi=yf[1];
        this.water=yf[2];
        this.gs=this.gsmax*this.water;
    }
    
    //runge-kutta
    public double[] rk4(double[] y0,double dt,Dfun dfun){
        double[] k1=dfun.getdf(y0);
        k1=scaleArray(k1,dt);
        double[] y1=addArray(y0,scaleArray(k1,0.5));
        double[] k2=dfun.getdf(y1);
        k2=scaleArray(k2,dt);
        double[] y2=addArray(y0,scaleArray(k2,0.5));
        double[] k3=dfun.getdf(y2);
        k3=scaleArray(k3,dt);
        double[] y3=addArray(y0,scaleArray(k3,0.5));
        double[] k4=dfun.getdf(y3);
        k4=scaleArray(k4,dt);
        double[] sumk=addArray(addArray(k1,k2),addArray(k3,k4));
        sumk=scaleArray(sumk,1/6);
        double[] yf=addArray(y0,sumk); 
        return yf;
    }
    
    public double[] scaleArray(double[] arr, double scaleFactor){
        double[] result;
        result=new double[arr.length];
        for (int i=0; i<arr.length; i++) {
            result[i] = arr[i] * scaleFactor;
        }
        return result;
    } 
    
    public double[] addArray(double[] arr1, double[] arr2){
        for (int i=0; i<arr1.length; i++) {
            arr1[i] = arr1[i] + arr2[i];
        }
        return arr1;
    }
    
    //get equilibrium signal 
    public double getSeq(){
        double smin = this.gsmin / this.gsmax;
        double seqval = ((1 + smin + this.alphag * this.irr)
                - sqrt(pow(1 + smin + this.alphag * this.irr,2)
                - 4 * this.thetag * (smin + this.alphag * this.irr))) / (2 * this.thetag);
        return seqval;
    }
    
    //get equilibrium vf, vj, vc
    private double getvxeq(double alphax, double thetax, double vmax, double vxmin){
        double fn = (((alphax * this.irr + 1 - vxmin)
                - sqrt(pow(alphax * this.irr + 1 - vxmin,2)  
                - 4 * alphax * this.irr * thetax * (1 - vxmin))) 
                / (2 * thetax)) + vxmin;
        double vxeq = vmax * fn;
        return vxeq;
    }
    
    private double getR(){
        double shit1 = (this.vj * (1 - this.poolT / this.tmax))
                * ((this.kt + this.poolT) / (this.vf * this.poolT));
        //no 5/3 in units of RuBP, corrected
        double pooR = (1 - shit1) * this.rmax;
        return pooR;       
    }
    
    private double geteqn2(double localci, double wc){
        double eqn2 = (((this.vf * this.poolT) / (this.kt + this.poolT))
                * (1 - this.poolR / this.rmax))
                - (1 + 2 * this.gamma / localci) * wc
                * (this.poolR / (this.kr + this.poolR));
        return eqn2;
    }
    
    private double getgg(){
        double gg = (this.gs * 1000000 * this.gb) / (this.gs * 1000000 + this.gb); 
        return gg;
    }
    
    private double getwc(double localci){
        double wc = (this.vc * pow(localci , 2)) / ((localci + this.ka)
                * (localci + this.kc * (1 + this.po2 / this.ko)));
        return wc;
    }
    
    private double getass1(double wc){
        double ass1 = wc * (this.poolR / (this.poolR + this.kr))
                - 0.5 * this.psi * this.poolG - this.rd;
        return ass1;
    }       
            
}
