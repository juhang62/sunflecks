/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package sunflecks;
import static java.lang.Math.*;

//comment which class from which package

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
    public double tauci = 79.8627652635799;
    public double taucd = 110.598886780825;
    public double taufi=471.475996185132;
    public double taufd=304.671929390908;
    
    
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

 

    public static void main(String[] args) { 
        Sunflecks obj=new Sunflecks(100);
        double[] time=new double[3600];
        double[] irr=new double[3600];
        for (int i=0;i<1800;i++){
            time[i]=i;
            irr[i]=500;
        }
        for (int i=1800;i<3600;i++){
            time[i]=i;
            irr[i]=100;
        }
        double[] ass=obj.rawrun(time,irr);
        System.out.println(ass[1]);  
    }
    
    public double[] rawrun(double[] time, double[] irr){
        //double dt=0.01;
        int nDatapts= time.length;
        double timeintv;
        double irrin;
        double[] assco2=new double[nDatapts];
        //Sunflecks simSunflecks=new Sunflecks(timeirr[0][1]); 
        for (int i=1;i<nDatapts;i++){
            timeintv=time[i]-time[i-1];
            irrin=(irr[i]+irr[i-1])/2;
            this.calcDyn(irrin,timeintv);
            //get An at each data point and write into output file
            assco2[i]=this.ass;
        }
        return assco2;
    }

    
    //constructor be called when instantiate the class
    public Sunflecks(double iniIRR) {
        this.calcSteady(iniIRR);
    }
    
    private String[] getOutput(Double curtime){
        String[] strArr= new String[15];
        strArr[0] =  String.valueOf(curtime);
        strArr[1]= String.valueOf(this.ass);
        strArr[2]= String.valueOf(this.irr);
        strArr[3]= String.valueOf(this.ci);
        strArr[4]= String.valueOf(this.vc);
        strArr[5]= String.valueOf(this.vf);
        strArr[6]= String.valueOf(this.vj);
        strArr[7]= String.valueOf(this.poolT);
        strArr[8]= String.valueOf(this.poolR);
        strArr[9]= String.valueOf(this.poolG);
        strArr[10]= String.valueOf(this.gs);
        strArr[11]= String.valueOf(this.signal);
        strArr[12]= String.valueOf(this.pi);
        strArr[13]= String.valueOf(this.water);
        strArr[14]= String.valueOf(this.seq);
        return strArr;
    }
    
    //private ici does not converge 315, 1242..
    public void calcSteady(double irr){
        this.irr=irr;
        this.vj=getvxeq(this.alphaj,this.thetaj,this.vjmax,0);
        this.vf=getvxeq(this.alphaf,this.thetaf,this.vfmax,this.vfmin);
        vc=getvxeq(this.alphac,this.thetac,this.vcmax,this.vcmin);
        updateStdygs();
        double newci=38.9; //first guess
        this.ci=39;
        int ici=0;
        while(abs(newci-this.ci) > 0.0001 && ici<1000) {
            ici = ici + 1;
            if (ici==999){
                System.out.println("ici error");
                System.out.println(irr);
            }
//            //double oldci=this.ci;
//            this.ci= newci;
//            double fx1=fci(this.ci);
//            double fx2=fci(this.ci-0.01);
//            //Netwon's Method: xnew=xold-f(xold)/f'(xold)
//            newci=this.ci-fx1*(0.01/(fx1-fx2)); 
            double oldci=this.ci;
            this.ci=newci;
            double fx1=fci(this.ci);
            double fx2=fci(oldci);
            //secent method
            newci=this.ci-fx1*((this.ci-oldci)/(fx1-fx2));
        }
        double gg = getgg(); 
        this.ass=(this.ca-this.ci)*gg/this.patm;                 
    }
    
    public void calcDyn(double irr, double timeinv){
        this.irr=irr;
        double dt=timeinv;
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
        //create a object of Differential Equation using my class Dfun
        Dfun photoDif=new Dfun(tau,this.rmax, this.tmax, this.kr,
                this.kt, this.gamma, this.psi,this.vj,vfeq,vceq, this.ci, wc);
        //initial value of state variables
        double[] y0={this.vc, this.vf, this.poolT, this.poolR, this.poolG};
        //calculate those values after one time step
        double[] yf=rkf45(y0,photoDif,dt); //use my rkf45 method 
        this.vc=yf[0];
        this.vf=yf[1];
        this.poolT=yf[2];
        this.poolR=yf[3];
        this.poolG=yf[4];
        wc=getwc(this.ci);
        this.ass=getass1(wc); 
        //System.out.println(this.ass);
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
            T2=this.tmax;
            //find t1
            for(int findt1_i=50; findt1_i<5000; findt1_i+=100){
                double Tportion=T2/(findt1_i+1);             
                for(int findt1_j=1; findt1_j<=findt1_i; findt1_j++){
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
            T1=this.tmax;
            //find t2
            for(int findt2_i=50; findt2_i<5000; findt2_i+=100){
                double Tportion=T1/(findt2_i+1);             
                for(int findt2_j=1; findt2_j<=findt2_i; findt2_j++){
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
    
    public double[] rkf45(double[] y0,Dfun dfun,double tspan){
        double dt=0.01;
        double dtmax=0.5;
        double dtmin=0.001;
        double t=0;
        double reltol=1e-4;
        double abstol=1e-4;
        int nvar=y0.length;        
        double[][] rkfinmd = new double[5][5];
        rkfinmd[0][0]=1.0/4;
        rkfinmd[1][0]=3.0/32;
        rkfinmd[1][1]=9.0/32;
        rkfinmd[2][0]=1932.0/2197;
        rkfinmd[2][1]=-7200.0/2197;
        rkfinmd[2][3]=7296.0/2197;
        rkfinmd[3][0]=439.0/216;
        rkfinmd[3][1]=-8;
        rkfinmd[3][2]=3680.0/513;
        rkfinmd[3][3]=-845.0/4104;
        rkfinmd[4][0]=-8.0/27;
        rkfinmd[4][1]=2;
        rkfinmd[4][2]=-3544.0/2565;
        rkfinmd[4][3]=1859.0/4104;
        rkfinmd[4][4]=-11.0/40;
        
        double[][] rkffi=new double[2][6];
        rkffi[0][0]=25.0/216;
        rkffi[0][1]=0;
        rkffi[0][2]=1408.0/2565;
        rkffi[0][3]=2197.0/4104;
        rkffi[0][4]=-1.0/5;
        rkffi[0][5]=0;
        rkffi[1][0]=16.0/135;
        rkffi[1][1]=0;
        rkffi[1][2]=6656.0/12825;
        rkffi[1][3]=28561.0/56430;
        rkffi[1][4]=-9.0/50;
        rkffi[1][5]=2.0/55;
    
        while (t < tspan) {
            if (t + dt > tspan){
                dt=tspan-t;
            } 
            double[][] k = new double[6][nvar];
            double[] y = y0;
            for (int i = 0; i < 6; i++) {
                k[i] = dfun.getdf(y);
                k[i]= scaleArray(k[i],dt);
                y = y0;
                if (i == 5) {
                    break;
                }
                for (int j = 0; j <= i; j++) {
                    y = addArray(y, scaleArray(k[j], rkfinmd[i][j]));
                }
            }

            double[] yf4 = y0;
            double[] yf5 = y0;
            for (int i = 0; i < 6; i++) {
                yf4 = addArray(yf4, scaleArray(k[i], rkffi[0][i]));
                yf5 = addArray(yf5, scaleArray(k[i], rkffi[1][i]));
            }
            double err = normmax(addArray(yf4,scaleArray(yf5,-1)));
            double ratio=err/(reltol*minabs(yf4)+abstol);
            if (ratio<1) {
                t = t + dt;
                y0 = yf4;
                dt=0.9*dt*pow(ratio,-0.2);
            }else{
                dt=0.9*dt*pow(ratio,-0.25);
            }
            
            if (dtmax < dt) {
                dt = dtmax;
            }
            if (dtmin > dt) {
                System.out.println("stiff");
            }
        }
        
        return y0;
    }
    
    public static double normmax(double[] vec){
        double max=0;
        for( int i = 0; i < vec.length; i++ ){
            if(abs(vec[i])>max){
                max=vec[i];
            }
        }
        return max;    
    }
    
    public static double minabs(double[] vec){
        double min=10000;
        for( int i = 0; i < vec.length; i++ ){
            if(abs(vec[i])<min){
                min=vec[i];
            }
        }
        return min;    
    }
    
    
    public static double[] scaleArray(double[] arr, double scaleFactor){
        double[] result;
        result=new double[arr.length];
        for (int i=0; i<arr.length; i++) {
            result[i] = arr[i] * scaleFactor;
        }
        return result;
    } 
    
    
    
    
    public static double[] addArray(double[] arr1, double[] arr2){
        double[] out=new double[5];
        for (int i=0; i<arr1.length; i++) {
            out[i] = arr1[i] + arr2[i];
        }
        return out;
    }
    
    public static double[] addMultiArray(double[] ... arr){
        int narr=arr.length;
        double[] sum=arr[1];
        for (int i=1; i<narr; i++) {
            sum = addArray(sum, arr[i]);
        }
        return sum;
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
