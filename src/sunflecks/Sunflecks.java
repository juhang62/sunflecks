/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package sunflecks;
import static java.lang.Math.*;
import au.com.bytecode.opencsv.CSVReader;
import au.com.bytecode.opencsv.CSVWriter;
import java.io.FileNotFoundException;

import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.StringWriter;
import java.util.List;

/**
 *
 * @author lfhan
 */
public final class Sunflecks {

    public double signal;
    public double pi;
    public double water;
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
    public double psi = 0.321430269302047;
    public double taugi = 215.5893;
    public double taugd = 324.2205;
    public double taupi = 206.6032;
    public double tauw = 208.7749;
    
    //psyn timeconst
    double tauci = 79.8627652635799;
    double taucd = 110.598886780825;
    double taufi=471.475996185132;
    double taufd=304.671929390908;
    
    
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
    public double time_interval;
    public double seq;

    
    

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws FileNotFoundException, IOException {
        Sunflecks simSunflecks=new Sunflecks(100);      
        double dt=0.01;
//        double[] timearr={
//            1, 3, 4, 5, 500};
//        double[] irrarr={
//        400, 500, 500,500,500};  
        
        double[] timearr;
        timearr=new double[563];
        double[] irrarr;
        irrarr=new double[563];
        CSVReader reader = new CSVReader(new FileReader("sim4488.csv"));
        CSVWriter writer = new CSVWriter(new FileWriter("output88.csv"), '\t');
        String [] nextLine;
        int iline=0;
        while ((nextLine = reader.readNext()) != null) {
            // nextLine[] is an array of values from the line
            timearr[iline]= Double.parseDouble(nextLine[0]);
            irrarr[iline]= Double.parseDouble(nextLine[2]);
            iline+=1;
        }
        
        double curtime=0;
        for (int i=0;i<timearr.length;i++){
            while (curtime<timearr[i]){
                curtime+=dt;
                simSunflecks.calcDyn(irrarr[i]);             
            }
            String str2 =  String.valueOf(simSunflecks.ass);
            String str1 =  String.valueOf(timearr[i]);
            String[] strarr= {str1,str2};
            writer.writeNext(strarr);
        }        
        writer.close();
        double teststr=simSunflecks.ass;
        System.out.println(curtime);
        System.out.println(teststr);
    }
    
    public Sunflecks(double iniIRR){
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
        double newci=this.ci+1;
        for(int ci_i=1; ci_i<500; ci_i++){
            double oldci=this.ci;
            this.ci=newci;
            double fx1=dynfci(this.ci);
            double fx2=dynfci(oldci);
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
        double[] tau={this.tauci,this.taucd,this.taufi,this.taufd};
        Dfun photoDif=new Dfun(tau,this.rmax, this.tmax, this.kr,
                this.kt, this.gamma, this.psi,this.vj,vfeq,vceq, this.ci, wc);
        double[] y0={this.vc, this.vf, this.poolT, this.poolR, this.poolG};
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
    
    public double dynfci(double localci){
        double wc=getwc(localci);
        double ass1 = getass1(wc);
        double gg = getgg();
        double ass2=(this.ca-localci)*gg/this.patm;
        return ass1-ass2;
    }
    
    public double fci(double localci){
        double wc = getwc(localci);
        this.poolT=this.tmax;
        this.poolR = getR(); 
        double eqn2 = geteqn2(localci,wc);       
        double T2 = 0;
        double T1 = 0;
        if (eqn2>0){
            T2=this.poolT;
            //find t1
            for(int findt1_i=1; findt1_i<500; findt1_i++){
                double Tportion=T2/pow(2,findt1_i);             
                for(int findt1_j=1; findt1_j<pow(2,findt1_i+1); findt1_j++){
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
                double Tportion=T1/pow(2,findt2_i);             
                for(int findt1_j=1; findt1_j<pow(2,findt2_i); findt1_j++){
                    T2=T1-Tportion*findt1_j;
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
    
    public void updateStdygs(){
        this.seq=getSeq();
        this.gs=this.seq*this.gsmax;
        this.signal=this.seq;
        this.pi=this.seq;
        this.water=this.seq;
    }
    
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
    
    public double getSeq(){
        double smin = this.gsmin / this.gsmax;
        double seqval = ((1 + smin + this.alphag * this.irr)
                - sqrt(pow(1 + smin + this.alphag * this.irr,2)
                - 4 * this.thetag * (smin + this.alphag * this.irr))) / (2 * this.thetag);
        return seqval;
    }
    
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
