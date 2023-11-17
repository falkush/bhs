package bhs;

import java.awt.AWTException;

import java.awt.image.BufferedImage;
import java.awt.image.DataBufferByte;
import java.io.File;
import java.io.IOException;

import javax.imageio.ImageIO;



public class bhs {
	static final int n=3;
	 public static void main(String[] args) throws AWTException, InterruptedException, IOException {
		
		final int height=1080;
		final int width=1920;
		
		int startm=0;
		int finm=1660;
		
		double schecker=0.1d;
		double roomsize=100d;
		double bhsize=100d*Math.PI*2d;
		
		double dist=1d;
		double sqsz=0.01d/6d;
		
		
		double[][] gmov = new double[finm][3];
		double[] gmov2 = new double[finm];
		double[] galpha = new double[finm];

		double[][][] gang = new double[finm][3][3];

		double[][] gspin = new double[finm][3];
		
		double lecossq;
		
		int i,j,pw,ph,l;
	
		double tsol;
		double tmin;
		int tmincoord;
		int setcoord;
		double[] coll = new double[2];
		
    	double d;
    	double tmpang;
    	double[] newpos = new double[3];
    	
    	double alpha;
    	double ttot;
    	double uhpy;
    	double dilf=1d;
    	
    	double tcont;
    	double sx,sy,sz;
    	double exitlgt;
    	
    	double[] xf = new double[3];
    	
    	double[] vecn2 = new double[3];
    	
    	double lerayon,lex1,lex2;
    	
    	
		int currentpix;
		double tf = 0;
	    
	    BufferedImage image = new BufferedImage(width,height,BufferedImage.TYPE_3BYTE_BGR);
	    
	    int nbframe=0;
	    
		double msqsz=-sqsz;
		double multy=((double)(1-width))*sqsz/2d;
		double multz=((double)(height-1))*sqsz/2d;
		
		double[] addy = new double[3];
		double[] addz = new double[3];
		double[] vectmp = new double[3];
		double[] vecn = new double[3];
		double[][] vecl = new double[height][width];
		
		int[] ctmp = new int[3];
			
		double[] pos = new double[3];
		double[] pos2 = new double[3];
		
		
		double exitangle2;
		
		int checker;
		double dc;
		
		double[][] x = new double[3][3];
		
		double dotp;
	
		double exitangle,lecos,leangle;
		
		
		double[][] newx = new double[3][3];
		
		double[] vec = new double[3];
	
			double discr,qa,qb,qc;
			
			double[] vecperp = new double[3];
				double vecperpn;
		
		
		xf[0]=-0.9999264065969468d;
		xf[1]=-0.011393387307123254d;
		xf[2]=-0.0041679870097623035d;
		

		for(i=0;i<3;i++) x[i][i]=1d;

		

		for(i=0;i<3;i++)
		{
			vec[i]=dist*x[0][i]+multy*x[1][i]+multz*x[2][i];
			addy[i]=sqsz*x[1][i];
			addz[i]=msqsz*x[2][i];
			pos[i]=0d;
		}

		
		for(ph=0;ph<height;ph++)
		{
			for(i=0;i<3;i++) vectmp[i]=vec[i];
			for(pw=0;pw<width;pw++)
			{
				vecl[ph][pw]=1d/veclgt(vec);
				for(i=0;i<3;i++) vec[i]+=addy[i];
			}
			for(i=0;i<3;i++) vec[i]=vectmp[i]+addz[i];
		}
	  
	
	    
	    final byte[] pixels =((DataBufferByte) image.getRaster().getDataBuffer()).getData();
	    
	    
	    gmov[0][0]=-10d;
	    
	    for(i=1;i<1200;i++)
		{
			gspin[i][1]=0.01d/4d;
			gspin[i][2]=0.003d/4d;
			gmov[i][0]=0.01d;
		}
	    
	    for(i=1200;i<1300;i++)
		{
			gmov[i][0]=0.01d;
		}
	    
	    for(i=1300;i<1600;i++)
	    {
	    	gang[i][0][1]=1/300d;
	    	gmov2[i]=0.01d;
	    }
	    
	    
	    
	    for(i=1600;i<1660;i++)
	    {
	    	gmov2[i]=0.01d;
	    	galpha[i]=(-500d/60d)*i+(1660d*500d/60d);
	    }
	   

       
	    while(nbframe<finm)
	    {
	    	
	    	for(i=0;i<3;i++)
	    	{
	    		 for(j=0;j<3;j++) pos[j]+=x[i][j]*gmov[nbframe][i]*dilf;
	    	}
	    	
	    	 for(j=0;j<3;j++) pos[j]+=xf[j]*gmov2[nbframe]*dilf;
	    	
	    	for(i=0;i<3;i++)
	    	{
	    		for(j=0;j<3;j++)
	    		{
	    			for(l=0;l<3;l++)
	    			{
	    				newx[i][l]=x[i][l]*Math.cos(gang[nbframe][i][j]*Math.PI)-Math.sin(gang[nbframe][i][j]*Math.PI)*x[j][l];
	    				newx[j][l]=x[j][l]*Math.cos(gang[nbframe][i][j]*Math.PI)+Math.sin(gang[nbframe][i][j]*Math.PI)*x[i][l];
	    				x[i][l]=newx[i][l];
	    				x[j][l]=newx[j][l];
	    			}
	    		}
	    	}

	    	
	    	////////////////////////////////////
	    	for(j=1;j<3;j++) {
	    		if(gspin[nbframe][j]!=0d) {
	    	d=veclgt(pos);
	    	  	

	    	
	    	tmpang=Math.PI-gspin[nbframe][j];
	    	
	    	for(i=0;i<3;i++)
	    	{
	    		newpos[i]=d*(Math.cos(tmpang)*x[0][i]+Math.sin(tmpang)*x[j][i]);
	    	}
	    	
	    	tmpang=gspin[nbframe][j];

	    	for(i=0;i<3;i++) pos[i]=newpos[i];
	    	
	    	for(l=0;l<3;l++)
			{
				newx[0][l]=x[0][l]*Math.cos(tmpang)-Math.sin(tmpang)*x[j][l];
				newx[j][l]=x[j][l]*Math.cos(tmpang)+Math.sin(tmpang)*x[0][l];
				x[0][l]=newx[0][l];
				x[j][l]=newx[j][l];
			}
	    	}}
	    	
	    	
		    ///////////////////////////////////////////		
	    	
	    	dilf=Math.sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]);
	    	if(dilf>1)dilf=1;
	    	
		    		if(nbframe>=startm) {
		    			
			    	currentpix=0;
			    	for(i=0;i<3;i++)
			    	{
			    		vec[i]=dist*x[0][i]+multy*x[1][i]+multz*x[2][i];
			    		addy[i]=sqsz*x[1][i];
			    		addz[i]=msqsz*x[2][i];
			    	}
			    	
			    	dc=Math.sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]);
				    	
			 	    for(ph=0;ph<height;ph++)
			 	    {
			 	    	for(i=0;i<3;i++) vectmp[i]=vec[i];
			 	    	for(pw=0;pw<width;pw++)
			 	    	{
		
			 	    		for(i=0;i<3;i++)
			 				{
			 					vecn[i]=vec[i]*vecl[ph][pw];
			 				}
			 	    		
			 	    		if(dc<1d)
			 	    		{
			 	    			dilf=dc;
			 	    			uhpy=1/dc;
			 	    			
			 	    			sx=pos[0]/dc;
			 	    			sy=pos[1]/dc;
			 	    			sz=pos[2]/dc;
			 	    			
			 	    			dotp=sx*vecn[0]+sy*vecn[1]+sz*vecn[2];
		 	    				leangle=Math.acos(dotp);
		 	    					
		 	    				vecperp[0]=vecn[0]-dotp*sx;
		 	    				vecperp[1]=vecn[1]-dotp*sy;
		 	    				vecperp[2]=vecn[2]-dotp*sz;
		 	    				
		 	    				vecperpn = Math.sqrt(vecperp[0]*vecperp[0]+vecperp[1]*vecperp[1]+vecperp[2]*vecperp[2]);
		 	    				
		 	    				vecperp[0]/=vecperpn;
		 	    				vecperp[1]/=vecperpn;
		 	    				vecperp[2]/=vecperpn;
		 	    				
		 	    				lecos=Math.cos(leangle);
		 	    				
		 	    				
		 	    				lecossq=lecos*lecos;
		 	    				
		 	    				
		 	    				
		 	    				lerayon=uhpy*Math.sqrt(1d/(1-lecossq));
		 	    				
		 	    				lex1=Math.sqrt(lerayon*lerayon - uhpy*uhpy)*Math.signum(lecos);
		 	    				lex2=Math.sqrt(lerayon*lerayon-1);
		 	    				
		 	    				
		 	    				exitlgt=lex2-lex1;
		 	    				exitangle=exitlgt%bhsize;
		 	    				exitangle=(2d*Math.PI/bhsize)*exitangle;
		 	    				
		 	    				exitangle2=Math.acos((1d/uhpy)*Math.sqrt(uhpy*uhpy-1d+lecossq));
		 	    				

		 	    				pos2[0]=Math.cos(exitangle)*sx+Math.sin(exitangle)*vecperp[0];
		 	    				pos2[1]=Math.cos(exitangle)*sy+Math.sin(exitangle)*vecperp[1];
		 	    				pos2[2]=Math.cos(exitangle)*sz+Math.sin(exitangle)*vecperp[2];

		 	    				
		 	    				vecn2[0]=Math.cos(exitangle+exitangle2)*sx+Math.sin(exitangle2+exitangle)*vecperp[0];
		 	    				vecn2[1]=Math.cos(exitangle+exitangle2)*sy+Math.sin(exitangle2+exitangle)*vecperp[1];
		 	    				vecn2[2]=Math.cos(exitangle+exitangle2)*sz+Math.sin(exitangle2+exitangle)*vecperp[2];
		 	    				
		 	    				tf=2d*Math.log((Math.sqrt((lex2-lex1)*(lex2-lex1)+(1-uhpy)*(1-uhpy))+Math.sqrt((lex2-lex1)*(lex2-lex1)+(1+uhpy)*(1+uhpy)))/(2d*Math.sqrt(uhpy)));
		 	    				tcont=0;
			 	    			
			 	    		}
			 	    		else
			 	    		{
			 	    			dilf=1d;
			 	    			qa=vecn[0]*vecn[0]+vecn[1]*vecn[1]+vecn[2]*vecn[2];
			 	    			qb=2d*(vecn[0]*pos[0]+vecn[1]*pos[1]+vecn[2]*pos[2]);
			 	    			qc=pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]-1d;
			 	    			
			 	    			discr=qb*qb-4d*qa*qc;
			 	    			
			 	    			if(discr<=0) {
				 	    			pos2[0]=pos[0];
				 	    			pos2[1]=pos[1];
				 	    			pos2[2]=pos[2];
				 	    			
				 	    			vecn2[0]=vecn[0];
				 	    			vecn2[1]=vecn[1];
				 	    			vecn2[2]=vecn[2];
				 	    			
				 	    			tcont=0;
				 	    			tf=0;
			 	    			}
			 	    			else
			 	    			{
			 	    				
			 	    				tcont=((-1d)*qb-Math.sqrt(discr))/(2d*qa);
			 	    				
			 	    				sx=vecn[0]*tcont+pos[0];
			 	    				sy=vecn[1]*tcont+pos[1];
			 	    				sz=vecn[2]*tcont+pos[2];
			 	    				
			 	    				dotp=sx*vecn[0]+sy*vecn[1]+sz*vecn[2];
			 	    				leangle=Math.acos(dotp);
			 	    					
			 	    				vecperp[0]=vecn[0]-dotp*sx;
			 	    				vecperp[1]=vecn[1]-dotp*sy;
			 	    				vecperp[2]=vecn[2]-dotp*sz;
			 	    				
			 	    				vecperpn = Math.sqrt(vecperp[0]*vecperp[0]+vecperp[1]*vecperp[1]+vecperp[2]*vecperp[2]);
			 	    				
			 	    				vecperp[0]/=vecperpn;
			 	    				vecperp[1]/=vecperpn;
			 	    				vecperp[2]/=vecperpn;
			 	    				
			 	    				lecos=Math.cos(leangle);
			 	    				lecos*=lecos;
			 	    			
			 	    				
			 	    				exitlgt=2d*(Math.sqrt(lecos/(1d-lecos)));
			 	    				exitangle=exitlgt%bhsize;
		
			 	    				exitangle=(2d*Math.PI/bhsize)*exitangle;
			 	    				

			 	    				pos2[0]=Math.cos(exitangle)*sx+Math.sin(exitangle)*vecperp[0];
			 	    				pos2[1]=Math.cos(exitangle)*sy+Math.sin(exitangle)*vecperp[1];
			 	    				pos2[2]=Math.cos(exitangle)*sz+Math.sin(exitangle)*vecperp[2];

			 	    				
			 	    				vecn2[0]=Math.cos(Math.PI+exitangle-leangle)*sx+Math.sin(Math.PI-leangle+exitangle)*vecperp[0];
			 	    				vecn2[1]=Math.cos(Math.PI+exitangle-leangle)*sy+Math.sin(Math.PI-leangle+exitangle)*vecperp[1];
			 	    				vecn2[2]=Math.cos(Math.PI+exitangle-leangle)*sz+Math.sin(Math.PI-leangle+exitangle)*vecperp[2];
			 	    				
			 	    		
			 	    				
			 	    				tf=2d*Math.log((exitlgt+Math.sqrt(exitlgt*exitlgt+4))/2d);
			 	    				
			 	    				
			 	    			}
			 	    		}
			 	    		
			 	    		
			 	    		tmin=(Math.signum(vecn2[0])*roomsize-pos2[0])/vecn2[0];
		 	    			tmincoord=0;
		 	    			
		 	    			tsol=(Math.signum(vecn2[1])*roomsize-pos2[1])/vecn2[1];
		 	    			if(tsol<tmin)
		 	    			{
		 	    				tmin=tsol;
		 	    				tmincoord=1;
		 	    			}
		 	    			
		 	    			tsol=(Math.signum(vecn2[2])*roomsize-pos2[2])/vecn2[2];
		 	    			if(tsol<tmin)
		 	    			{
		 	    				tmin=tsol;
		 	    				tmincoord=2;
		 	    			}
		 	    			
		 	    			setcoord=0;
		 	    			for(i=0;i<3;i++)
		 	    			{
		 	    				if(i!=tmincoord)
		 	    				{
		 	    					coll[setcoord]=pos2[i]+tmin*vecn2[i];
		 	    					setcoord++;
		 	    				}
		 	    			}
		 	    			
		 	    			checker=((int)Math.floor(coll[0]*schecker))%2;
		 	    			checker+=((int)Math.floor(coll[1]*schecker))%2;
		 	    			if(checker<0) checker+=2;
		 	    			checker%=2;
		 	    			
		 	    			
		 	    			if(tmincoord==2) {
		 	    			ctmp[0]=255;
		 	    			ctmp[1]=255;
		 	    			ctmp[2]=255; }
		 	    			else if(tmincoord==0)
		 	    			{
		 	    				if(vecn2[0]<0) ctmp=rnbw((1d/8d)-coll[0]/(8d*roomsize));
		 	    				else ctmp=rnbw((1d/2d)+(1d/8d)+coll[0]/(8d*roomsize)); 
		 	    			}
		 	    			else
		 	    			{
		 	    				if(vecn2[1]<0) ctmp=rnbw((1d/4d)+(1d/8d)+coll[0]/(8d*roomsize));
		 	    				else ctmp=rnbw((3d/4d)+(1d/8d)-coll[0]/(8d*roomsize));
		 	    			}
		 	    			
		 	    			if(galpha[nbframe]!=0)
		 	    			{
		 	    			  ttot=tf+tcont+tmin;
		 	    			  alpha=((-1d)/galpha[nbframe])*ttot+1;
		 	    			  if(alpha<0) alpha=0d;
		 	    			}
		 	    			else alpha=1d;
		 	    			
		 	    			if(checker==0) {
		 	    		
				 	    		pixels[currentpix]=(byte)0;
				 	    		currentpix++;
				 	    		pixels[currentpix]=(byte)0;
				 	    		currentpix++;
				 	    		pixels[currentpix]=(byte)0;
				 	    		currentpix++;
		 	    			}
		 	    			else
		 	    			{
		 	    				pixels[currentpix]=(byte)(alpha*ctmp[0]);
				 	    		currentpix++;
				 	    		pixels[currentpix]=(byte)(alpha*ctmp[1]);
				 	    		currentpix++;
				 	    		pixels[currentpix]=(byte)(alpha*ctmp[2]);
				 	    		currentpix++;
		 	    			}
		 	    		
			 
			 	    		
			 	    		for(i=0;i<3;i++) vec[i]+=addy[i];
			 	    	}
			 	    	
			 	    	for(i=0;i<3;i++) vec[i]=vectmp[i]+addz[i];
			 	    }
			 	    	
			 	    	ImageIO.write(image, "bmp", new File("C:/test/bhs/img"+String.format("%04d", nbframe)+".bmp/"));
					      System.out.println(nbframe);

			 	    }
			 	    			
	    	nbframe++;
	    }
	}
	 
	 public static double veclgt(double[] v)
		{
			return Math.sqrt(vecprod(v,v));
		}
	 
	 public static double vecprod(double[] v1, double[] v2)
		{
			double ret=0;
			for(int i=0;i<n;i++) ret+=v1[i]*v2[i];
			return ret;
		}

	
	public static int[] rnbw(double x)
	{
		int[] ret = new int[3];
		double tmp;
		
		tmp=x%(1d/6d);
		
		if(x<1d/6d)
		{
			ret[0]=255;
			ret[1]=(int) (1530d*tmp);
		}
		else if(x<1d/3d)
		{
			ret[1]=255;
			ret[0]=(int) (255d-1530d*tmp);
		}
		else if(x<0.5d)
		{
			ret[1]=255;
			ret[2]=(int) (1530d*tmp);
		}
		else if(x<2d/3d)
		{
			ret[2]=255;
			ret[1]=(int) (255d-1530d*tmp);
		}
		else if(x<5d/6d)
		{
			ret[2]=255;
			ret[0]=(int) (1530d*tmp);
		}
		else
		{
			ret[0]=255;
			ret[2]=(int) (255d-1530d*tmp);
		}
		
		return ret;
	}
	

		

}