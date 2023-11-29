public class Delphi {


	//qlog.h
	boolean iautocon,iper[]=new boolean[3],iconc,ibios,isite,iatout,diff,isph;
	boolean ipdbwrt,ifrcwrt,ifrcrd,ipoten,igraph,imem,ihome,icheb;
	public boolean phiwrt,logs,logc,loga,logg,ipdbrd,logas,epswrt,iacent;
	boolean isch,isen,iacs,irea,ibufz,inrgwrt,iexun,iwgcrg,iuspec;
	boolean idbwrt,ixphiwt,iyphiwt,izphiwt,ixphird,iyphird,izphird;
	boolean isita,isitq,isitp,isitf,isitr,isitc,isitx,isiti,iself;
	boolean isitrf,isitcf,isitt;
	public boolean isolv,isrf,ibem;
 
	String epsnam,phinam,frcnam,mpdbnam,updbnam,ufrcnam,srfnam; 
	String centnam,siznam,crgnam,pdbnam,frcinam,phiinam;
	String prmnam,usernam,uhomenam,scrgnam,nrgnam,gcrgnam;
	String dbnam,xphonam,yphonam,zphonam,xphinam,yphinam,zphinam;
	String toplbl;
 
	float scale,prbrad,exrad,perfil,rionst,repsout,repsin,gten,uspec;
	float offset[]=new float[3],epsout,epsin,radprb,acent[]=new float[3],epkt,fpi,deblen;
	int igrid,nlit,nnit,ibctyp,icon1,icon2;
	int bufz[]=new int[6];
	int epslen,philen,srflen,frclen,mpdblen,updblen,ufrclen;
	int phifrm,epsfrm,frcfrm,mpdbfrm,updbfrm,ufrcfrm;
	int centlen,sizlen,crglen,pdblen,frcilen,phiilen;
	int centfrm,sizfrm,crgfrm,pdbfrm,frcifrm,phiifrm;
	int prmfrm,scrgfrm,scrglen,prmlen,userlen,uhomelen;
	int nrgfrm,nrglen,gcrglen,gcrgfrm,dblen,xphilen,yphilen;
	int zphilen,xphipos,yphipos,zphipos,xphopos,yphopos,zphopos;
	int xpholen,ypholen,zpholen;

  	//acc2.h
	
	int idmax=50;
	float r0[],r02[],rs2[];
	int lcb,mcb,ncb;
	float xo,yo,zo,cbai;
	int lcb1,mcb1,ncb1;
	float grdi,mnx,mny,mnz;
	int ast[],extot;
	float expos[],tary[]=new float[2];
	float cmin[]=new float[3];
	float cmax[]=new float[3];
	float rdmx;

	//pointer.h
	int iepsmp[],idebmap[],ioff[],idpos[];
	float phimap[],phimap1[];
	float phimap2[],phimap3[];
	float db[],sf1[],sf2[];
	float qmap1[],qmap2[];
	float debmap1[],debmap2[];
	float bndx1[],bndx2[],bndx3[],bndx4[];
	float bndx[],bndy[],bndz[];
	int   neps[],keps[];
	float cgrid[],spdiv[],sen[];
	float spot[],sqs[];
	int   iepsv[];
	float rfield[];
	int   pls[];
	float ast2[];
	int   ibnd[];
	int   ibgrd[];
      	int   bndeps[];
	int   cbn1[],cbn2[],cbal[];
	int   icbn[];
 	int   atndx[],iab1[],iab2[],icume[];
	int   iexpos[],atsurf[];
	float scspos[];
	float atpos[],xn2[],rad3[];
	float chrgv4[],atinf[];
	float qval[],gchrg[];
	int   gchrgp[];
	int iqpos[];
	float chrgv2[],cqs[];
	float vert[];
	int   vtemp[];
	float scsnor[],vnorm[];
	int   nsel[],vindx[];
	int vtlen[],vtlst[];
	int tmlst[],vtpnt[];
	float schrg[];
	int   crgatn[];
	float phimap4[];
	 


	//qdiffpar5.h
	int nclist = 1500;
	int nrlist = 1500;
	int nrmax = 1500; 
	int ncmax = 1500 ; 
	int ncrgmx = 50000 ; 
	int ngrid = 257;  
	int ngp = 274626;
	int nhgp = 137313;
	int natmax = 50000; 
	int ngcrg = 80000;
	int nbgp = 4226;
	int nsp = 40000;

	int limeps[]=new int[6];
	float atmcrg[],chgpos[],cgbp[],gval[];
	float oldmid[]=new float[3],oldmid1[]=new float[3];
	boolean logtab[]=new boolean[100];
	float scale1,rmmin,rmmax;
	int ibc;

	
		
 	//my declaration 
	long finish=0;
	int natom=0;
	float cran[]=new float[3];
	Molecule molecule=null;
	float rmaxdim=0;
	float cqplus[]=new float[3];
	float cqmin[]=new float[3];
	float qmin,qnet,qplus;
	int   ibnum;
	int   nqass,nqgrd;
	int   vtot,itot	;
	int   icount2a,icount2b,icount1a,icount1b;
	float om1,om2;
	public Surface surface;
	public Delphi(Molecule m) {
		 defprm();
		 molecule=m;
	}
 
	public static long cputime(long start) {

		return  System.currentTimeMillis()-start; 
	}
	
	String getTimeStr() {
		Date s=new Date ();
		return s.toString();	
	}

	void wrt(int ipgn){
	 

		if(ipgn==0) return;
		System.out.println("   ");
		System.out.println(" ___________________DelPhi II____________________  ");
		System.out.println("|                                                | ");
		System.out.println("| A program to solve the PB equation             | ");
		System.out.println("| in 3D, using non-linear form, incorporating    | ");
		System.out.println("| 2 dielectric regions, ionic strength, periodic | ");
		System.out.println("| and focussing boundary conditions, utilizing   | ");
		System.out.println("| stripped optimum successive over-relaxation    | ");
		System.out.println("| and an improved algorithm for mapping the      | ");
		System.out.println("| Mol. Surface to the finite-Difference grid     | ");
		System.out.println("|__________________          ____________________| ");
		System.out.println("                    DelPhi II                      ");
		System.out.println("  ");		 
		String day= getTimeStr();
        	System.out.println("program started on: "+day);          
		 
	}
	

	void defprm() {
 
		// boolean parameters used in qdiff, mainly from parameter file
	
	 
 
		//     data conversion to kT/charge at 25 celsius
 
     
		epkt=561.0f;
		fpi=12.566f;
 
		toplbl="qdiffxas: qdiffxs4 with an improved surfacing routine";
		isch=false;
		isen=false;
		iacs=false;
		irea=false;
		iautocon=true;
		iper[0]=false;
		iper[1]=false;
		iper[2]=false;
		iconc=false;
		ibios=false;
		isite=false;
		iatout=false;
		diff=false;
		isph=false;
		phiwrt=false;
		logs=false;
		logc=false;
		loga=false;
		logg=false;
		logas=false;
		ipdbrd=false;
		epswrt=false;
		iacent=false;
		ipdbwrt=false;
		ifrcwrt=false;
		ifrcrd=false;
		ipoten=false;
		igraph=false;
		imem=false;
		ihome=false;
		ibufz=false;
		inrgwrt=false;
		iwgcrg=false;
		iuspec=false;
		icheb=false;
		idbwrt=false;
		ixphird=false;
		iyphird=false;
		izphird=false;
		ixphiwt=false;
		iyphiwt=false;
		izphiwt=false;
		isita=false;
		isitq=false;
		isitp=false;
		isitf=false;
		isitr=false;
		isitc=false;
		isitx=false;
		isiti=false;
		iself=false;
		isitrf=false;
		isitcf=false;
		isitt=false;
		isolv=true;
		isrf=false;
		//
		icon1=10;
		icon2=1;
		epsnam="fort.17";
		phinam="fort.14";
		srfnam="grasp.srf";
		frcnam="fort.16";
		mpdbnam="fort.19";
		updbnam="fort.20";
		ufrcnam="fort.21";
		centnam="fort.15";
		pdbnam="fort.13";
		crgnam="fort.12";
		siznam="fort.11";
		phiinam="fort.18";
		frcinam="fort.15";
		prmnam="fort.10";
		scrgnam="scrg.dat";
		nrgnam="energy.dat";
		gcrgnam="crg.dat";
		dbnam="db.dat";
		xphinam="fort.31";
		yphinam="fort.32";
		zphinam="fort.33";
		xphonam="fort.34";
		yphonam="fort.35";
		zphonam="fort.36";
		//
		epslen=7;
		philen=7;
		srflen=9;
		frclen=7;
		mpdblen=7;
		updblen=7;
		ufrclen=7;
		centlen=7;
		pdblen=7;
		crglen=7;
		sizlen=7;
		phiilen=7;
		frcilen=7;
		scrglen=8;
		nrglen=10;
		prmlen=7;
		gcrglen=7;
		dblen=6;
		xphilen=7;
		yphilen=7;
		zphilen=7;
		xpholen=7;
		ypholen=7;
		zpholen=7;
		//
		phifrm=0;
		epsfrm=0;
		frcfrm=0;
		mpdbfrm=0;
		updbfrm=0;
		ufrcfrm=0;
		pdbfrm=0;
		crgfrm=0;
		sizfrm=0;
		prmfrm=0;
		phiifrm=0;
		frcifrm=0;
		scrgfrm=0;
		nrgfrm=0;
		gcrgfrm=0;
		scale=10000.f;
		radprb=1.4f;
		exrad=0.0f;
		perfil=10000.f;
		rionst=0.0f;
		repsout=80f;
		epsout=80f;
		repsin=2f;
		epsin=2f;
		gten=0.0f;
		uspec=0.9975f;
		offset[0]=0.f;
		offset[1]=0.f;
		offset[2]=0.f;
		acent[0]=0.0f;
		acent[1]=0.0f;
		acent[2]=0.0f;
		int i;
		for( int j=1;j<=3;j++){
			i=newIndexTwo(1,j,3);
			bufz[i]=0;
			i=newIndexTwo(2,j,3);
			bufz[i]=0;
		}
		igrid=0;
		nlit=0;
		nnit=0;
		ibctyp=2;
		//
	}

	final static int mod(int a,int b) {
	
		return (a-a/b*b);

        }

	final static int newIndexOne( int a) {
		return a-1;
	}

	final static int newIndexOne( int a,int b) {
		return a-b;
	}
	
	final static int newIndexTwo(int a,int b,int nb) {
		return (a-1)*nb+b-1;
	}

	final static int newIndexTwo(int a,int b,int na,int nb,int blen) {
		return (a-na)*blen+b-nb;
	}
	
	final static int newIndexThree(int a,int b,int c,int n) {
		
		return  (a-1)*n*n+(b-1)*n+c-1;
	}
	
	final static int newIndexFour(int a,int b,int c,int d,int n) {
		
		return  3*(a-1)*n*n+3*(b-1)*n+3*(c-1)+d-1;
	}

	final static int newIndexThree(int a,int b,int c,int na,int n) {
		
		return  (a-na)*n*n+(b-na)*n+c-na; 
	}


  	public void setPeriodic(boolean x,boolean y,boolean z) {
		iper[0]=x;
		iper[1]=y;
		iper[2]=z;
	}
	public void setBoundary(int a) {
		ibctyp=a;
	}
	 
	public void setPrm() {
		setPrm(65,70.0f, 2.f,80.f,1.6f,0f,0f);
	}

	public void setPrm(int sc,float per,float rin,float rout,float prb,float exr,float ron) {
		igrid=sc;
		perfil=per;
		repsin=rin;
		repsout=rout;
		radprb=prb;
		exrad=exr;	
		rionst=ron;
		float dfact = 3.047f;
		if((repsin<0.)||(repsout<0.)){
       			repsin=Math.abs(repsin);
        		repsout=Math.abs(repsout);
        		diff=true;
        	}
		if(rionst>0.) {
          		deblen = (float)(dfact/Math.sqrt(rionst));
        	}else {
          		deblen = 1.e6f;
       		}
		epsin = repsin/epkt;
        	epsout = repsout/epkt;
	}
	public void setPrm(float sc,float per,float rin,float rout,float prb,float exr,float ron) {
	 
		scale=sc;
		perfil=per;
		repsin=rin;
		repsout=rout;
		radprb=prb;
		exrad=exr;	
		rionst=ron;
		float dfact = 3.047f;
		if((repsin<0.)||(repsout<0.)){
       			repsin=Math.abs(repsin);
        		repsout=Math.abs(repsout);
        		diff=true;
        	}
		if(rionst>0.) {
          		deblen = (float) (dfact/Math.sqrt(rionst));
        	}else {
          		deblen = 1.e6f;
       		}
		epsin = repsin/epkt;
        	epsout = repsout/epkt;

	}
  

	void getatm() {
 
                
		int n=0;
		int i;
		natom=0;
                for(Molecule m=molecule;m!=null;m=m.next)
                for(Strand s=m.strand;s!=null;s=s.next)
                for(Residue r=s.residue;r!=null;r=r.next)
                for(Atom a=r.atom;a!=null;a=a.next) {
                        natom++;
                }
		
                if(atpos==null)  atpos=new float[natom*3];
                if(rad3==null)   rad3=new float[natom];
                if(chrgv4==null) chrgv4=new float[natom];
                
                for(Molecule m=molecule;m!=null;m=m.next)
                for(Strand s=m.strand;s!=null;s=s.next)
                for(Residue r=s.residue;r!=null;r=r.next)
                for(Atom a=r.atom;a!=null;a=a.next) {
                        n++;
			i=newIndexTwo(n,1,3);
                        atpos[i]= a.coord[0];
                        atpos[i+1]= a.coord[1];
                        atpos[i+2]= a.coord[2];
                        rad3[n-1]=a.radius;
                        chrgv4[n-1]=a.charge;
                }
		
        }

 	
	

	void ex(int egrid[], int ivn,int iv,int ivert[],float vc[],String vdat1,int lvdat1){
 
        	int mxtri=200000;
       		int mxvtx=100000; 
		int e2[]=new int[257*257*2];
		int ipv[]=new int[5*257*257];					
		int ivslb[]=new int[3*32000];
      		String uplbl;
      		String nxtlbl;
      		String toplbl;
      		String botlbl;
      		float x2,y2,z2;
      		int je,ix,iy,iz,ibox,j,n,k,i,i1,indx,i2;	
 		int itind[]=DelphiConst.itind, itrn[]=DelphiConst.itrn;
		
 
		float fourth=1.0f/4.0f;
	 	
 
		//**open[11,file="/usr/local/bin/v3.dat", form="unformatted");
		//**read[11]itind,itrn;
		//**close[11]; 		 
		// loop for times
 
		for( iy=1;iy<=igrid;iy++) {
	    		for( ix=1;ix<=igrid;ix++) {
				je=(ix-1)*257*5+(iy-1)*5+1-1;
				ipv[je]=0;
				ipv[je+1]=0;
				ipv[je+2]=0;
				ipv[je+3]=0;
				ipv[je+4]=0;
	  		}
		}
 
		for( iy=1;iy<=igrid;iy++) {
	    		for( ix=1;ix<=igrid;ix++) {
				je=newIndexFour(ix,iy,1,1,igrid);
				//System.out.println("iy,ix "+iy+" "+ix+" "+(egrid[je]));
				if(egrid[je]<=0) {
		    			je=(ix-1)*257*2+(iy-1)*2+2-1; 
		    			e2[je]=0;
				}else {
					je=(ix-1)*257*2+(iy-1)*2+2-1; 	
		    			e2[je]=1;
				}
	    		}
		}
		int iii=0;
 		for( iy=1;iy<=igrid;iy++) {
	    		for( ix=1;ix<=igrid;ix++) {
				for( iz=1;iz<=igrid;iz++) { 
					je=newIndexFour(ix,iy,iz,1,igrid);
					iii+=egrid[je];
		}}}
		//System.out.println("iii "+iii);
		
		ibox=0;
		iv=0;
		ivn=0;
		for( iz=1;iz<=igrid-1;iz++) { 
			for( iy=1;iy<=igrid;iy++) {
	   		for( ix=1;ix<=igrid;ix++) {
				je=(ix-1)*257*2+(iy-1)*2+1-1; 
				e2[je]=e2[je+1];
				je=newIndexFour(ix,iy,iz+1,1,igrid);
				if(egrid[je]<=0) {
					je=(ix-1)*257*2+(iy-1)*2+2-1; 
					e2[je]=0;
				}else {
					je=(ix-1)*257*2+(iy-1)*2+2-1; 
					e2[je]=1;
				}
	    		}
			}
 
			k=0;
			for( iy=1;iy<=igrid-1;iy++) {
	 
				j=e2[(iy-1)*2]+e2[iy*2]+e2[iy*2+1]+e2[(iy-1)*2+1];
				//System.out.println("iy,iz"+" "+iy+" "+iz+" aa "+j);
				for( ix=1;ix<=igrid-1;ix++) {
 
					i2=e2[ix*2*257+(iy-1)*2] + 
						e2[ix*2*257+iy*2]+ e2[ix*2*257+iy*2+1 ] + 
						e2[ix*2*257+(iy-1)*2+1];
					i=i2+j;
					j=i2;
 
					if((i!=0)&&(i!=8)) {
	    					k=k+1;
 
						// determine index of box, 1,254
						// it was found to be better to calculate the index NOW, not later
 
	    					indx=e2[(ix-1)*2*257+(iy-1)*2]+
							2*e2[(ix-1)*2*257+(iy)*2 ] + 
							4*e2[(ix)*2*257+(iy)*2]
							+8*e2[(ix)*2*257+(iy-1)*2  ]+ 
							16*e2[(ix-1)*2*257+(iy-1)*2+1]
							+32*e2[(ix-1)*2*257+(iy)*2+1]
							+ 64*e2[(ix)*2*257+(iy)*2+1]
							+128*e2[(ix)*2*257+(iy-1)*2+1];
							je=newIndexTwo(k,1,3);
							ivslb[je]=ix;
							ivslb[je+1]=iy;
							ivslb[je+2]=indx;							
					}
					// loop over those verticies in the triangle list for this index
 
					// loop to next box
				}
			}
 
			z2=(float)(2*iz);
			for( i=1;i<=k;i++) {
				je=newIndexTwo(i,1,3);
				ix=ivslb[je];
				iy=ivslb[je+1];
				indx=ivslb[je+2];
 				//System.out.println("iw...."+i+" "+ix+" "+iy+" "+indx+" "+k);
				x2=(float)(2*ix);
				y2=(float)(2*iy);
				for(i1=itind[indx];i1<=itind[indx+1]-1;i1++) {
					iv=iv+1;
					//System.out.println("i1..."+i1+" "+(itind[indx])+" "+indx);
					n=itrn[i1-1];
 
					// for each vertex, if it has a number in the vertex slab use it
					// if not, increment the new vertex counter anf fill ivc with the
					// right coordinaates
 
					//goto[10,20,30,40,50,60,70,80,90,100,110,120] n;
					
					if(n<=1||n>12) {
						je=(ix-1)*257*5+(iy-1)*5+3-1;
				 		j=ipv[je];
						if(j!=0) {
							ivert[iv-1]=j;
						}else {
	    						ivn=ivn+1;
							
	    						ipv[je]=ivn;
	    						ivert[iv-1]=ivn;
							je=newIndexTwo(ivn,1,3);
	    						vc[je]=x2;
	    						vc[je+1]=y2;
	    						vc[je+2]=z2+1;
						}
						//System.out.println("nnnn "+n+" "+j+" "+i1+" "+(itind[indx]));
					}
 					else if(n==2) {
						je=(ix-1)*257*5+(iy-1)*5+5-1;
						j=ipv[je];
						if(j!=0) {
	    						ivert[iv-1]=j;
						}else {
	    						ivn=ivn+1;
							
	    						ipv[je]=ivn;
	    						ivert[iv-1]=ivn;
							je=newIndexTwo(ivn,1,3);
	    						vc[je]=x2;
	    						vc[je+1]=y2+1;
	    						vc[je+2]=z2+2;
						}
					}	
					else if(n==3) {
						je=(ix-1)*257*5+(iy+1-1)*5+3-1;
 						j=ipv[je];
						if(j!=0) {
	    						ivert[iv-1]=j;
						}else {
	    						ivn=ivn+1;
							
	    						ipv[je]=ivn;
	    						ivert[iv-1]=ivn;
							je=newIndexTwo(ivn,1,3);
	    						vc[je]=x2;
	    						vc[je+1]=y2+2;
	    						vc[je+2]=z2+1;
						}
					}
 
					else if(n==4) {
						je=(ix-1)*257*5+(iy-1)*5+2-1;
						j=ipv[je];
						if(j!=0) {
							ivert[iv-1]=j;
						}else {
							ivn=ivn+1;
							
							ipv[je]=ivn;
							ivert[iv-1]=ivn;
							je=newIndexTwo(ivn,1,3);
							vc[je]=x2;
							vc[je+1]=y2+1;
							vc[je+2]=z2;
						}
					}
					else if(n==5) {
						je=(ix-1)*257*5+(iy-1)*5+1-1;
 						j=ipv[je];
						if(j!=0) {
							ivert[iv-1]=j;
						}else {
							ivn=ivn+1;
							
							ipv[je]=ivn;
							ivert[iv-1]=ivn;
							je=newIndexTwo(ivn,1,3);
							vc[je]=x2+1;
							vc[je+1]=y2;
							vc[je+2]=z2;
						}
					}
					else if(n==6) {
 						je=(ix-1)*257*5+(iy-1)*5+4-1;
 						j=ipv[je];
						if(j!=0) {
							ivert[iv-1]=j;
						}else {
							ivn=ivn+1;
							
							ipv[je]=ivn;
							je=newIndexTwo(ivn,1,3);
							ivert[iv-1]=ivn;
							vc[je]=x2+1;
							vc[je+1]=y2;
							vc[je+2]=z2+2;
						}
					}
	 				else if(n==7) {
 						je=(ix-1)*257*5+(iy+1-1)*5+4-1;
 						j=ipv[je];
						if(j!=0) {
							ivert[iv-1]=j;
						}else {
							ivn=ivn+1;
							ipv[je]=ivn;
							ivert[iv-1]=ivn;
							je=newIndexTwo(ivn,1,3);
							vc[je]=x2+1;
							vc[je+1]=y2+2;
							vc[je+2]=z2+2;
						}
					}	
					else if(n==8) {
 						je=(ix-1)*257*5+(iy+1-1)*5+1-1;
						j=ipv[je];
						if(j!=0) {
							ivert[iv-1]=j;
						}else {
							ivn=ivn+1;
							ipv[je]=ivn;
							je=newIndexTwo(ivn,1,3);
							ivert[iv-1]=ivn;
							vc[je]=x2+1;
							vc[je+1]=y2+2;
							vc[je+2]=z2;
						}
					}
					else if(n==9) {
 						je=(ix+1-1)*257*5+(iy-1)*5+3-1;
 						j=ipv[je];
						if(j!=0) {
							ivert[iv-1]=j;
						}else {
							ivn=ivn+1;
							ipv[je]=ivn;
							ivert[iv-1]=ivn;
							//System.out.println("ivn.."+ivn+" "+ix+" "+iy+" "+i1);
							je=newIndexTwo(ivn,1,3);
							vc[je]=x2+2;
							vc[je+1]=y2;
							vc[je+2]=z2+1;
						}
					}	
					else if(n==10) {
 						je=(ix+1-1)*257*5+(iy-1)*5+5-1;
 						j=ipv[je];
						if(j!=0) {
							ivert[iv-1]=j;
						}else {
							ivn=ivn+1;
							ipv[je]=ivn;
							ivert[iv-1]=ivn;
							je=newIndexTwo(ivn,1,3);
							vc[je]=x2+2;
							vc[je+1]=y2+1;
							vc[je+2]=z2+2;
					}
					}
					else if(n==11) {
						je=(ix+1-1)*257*5+(iy+1-1)*5+3-1;
						j=ipv[je];
						if(j!=0) {
							ivert[iv-1]=j;
						}else {
							ivn=ivn+1;
							ipv[je]=ivn;
							ivert[iv-1]=ivn;
							je=newIndexTwo(ivn,1,3);
							vc[je]=x2+2;
							vc[je+1]=y2+2;
							vc[je+2]=z2+1;
						}
					}
					else if(n==12) {
 						je=(ix+1-1)*257*5+(iy-1)*5+2-1;
 						j=ipv[je];
						if(j!=0) {
							ivert[iv-1]=j;
						}else {
							ivn=ivn+1;
							ipv[je]=ivn;
							ivert[iv-1]=ivn;
							je=newIndexTwo(ivn,1,3);
							vc[je]=x2+2;
							vc[je+1]=y2+1;
							vc[je+2]=z2;
						}
					}
	
				}
			}
 
			for( i=1;i<=k;i++) {
				je=newIndexTwo(i,1,3);
				ix=ivslb[je];
				iy=ivslb[je+1];
				je=(ix-1)*257*5+(iy-1)*5+1-1;
				ipv[je+2]=0;
				ipv[je]=ipv[je+3];
				ipv[je+1]=ipv[je+4];
				ipv[je+3]=0;
				ipv[je+4]=0;
			}
 
		}
 
		// it loop for timing
		finish= cputime(finish);
		System.out.println("number of vertices= "+ivn);
		System.out.println("number of triangles = "+iv/3);
		vtot=ivn;
		itot=iv;
		
	}
 	

	void extrm() {
 
		// find extrema and calculate scale according to them and
		// to the percent box fill
 		int j=0;
		cmin[0]=6000;
		cmin[1]=6000;
		cmin[2]=6000;
		cmax[0]=-6000;
		cmax[1]=-6000;
		cmax[2]=-6000;
		
		for( int ix=1;ix<=natom;ix++) {
			j=newIndexTwo(ix,1,3);
			cmin[0]=Math.min(cmin[0],atpos[j]-rad3[ix-1]);
			cmin[1]=Math.min(cmin[1],atpos[j+1]-rad3[ix-1]);
			cmin[2]=Math.min(cmin[2],atpos[j+2]-rad3[ix-1]);
			cmax[0]=Math.max(cmax[0],atpos[j]+rad3[ix-1]);
			cmax[1]=Math.max(cmax[1],atpos[j+1]+rad3[ix-1]);
			cmax[2]=Math.max(cmax[2],atpos[j+2]+rad3[ix-1]);
		} 
		oldmid[0]=(cmax[0]+cmin[0])/2;
		oldmid[1]=(cmax[1]+cmin[1])/2;
		oldmid[2]=(cmax[2]+cmin[2])/2;
		cran[0]=cmax[0]-cmin[0];
		cran[1]=cmax[1]-cmin[1];
		cran[2]=cmax[2]-cmin[2];
		rmaxdim=Math.max(cran[0],cran[1]);
		rmaxdim=Math.max(rmaxdim,cran[2]);
	 	
	}
	
	void wrtprm() {

		String bclab[]= { "zero","dipolar", "focussing","coulombic"};

		System.out.println(" ");

		System.out.println("grid size:"+igrid);
		if(scale<1.e-6) {
			System.out.println("percent of box to be filled:"+perfil);
		}
		else {
			System.out.println("scale,in grids/A, set to be:"+scale);
		}
		System.out.println("object centred at (gu):"+offset[0]+" "+offset[1]+" "+offset[2]);
		if(isolv){
			System.out.println("inner,outer dielectrics    :"+repsin+" "+repsout);
			System.out.println("ionic strength (M)         :"+rionst);
			System.out.println("debye length (A)           :"+deblen);
			System.out.println("ion exclusion radius (A)   :"+exrad);
			System.out.println("probe radius (A)           :"+radprb);
			System.out.println("boundary conditions        :"+bclab[ibctyp]);
			System.out.println("x,y,z periodic bc. flags   :"+iper[0]+" "+iper[1]+" "+iper[2]);
			if(iautocon) {
				if(gten>0.){
					System.out.println("convergence by grid energy:"+gten+" kt");
				}
				else {
					System.out.println("# of linear iterations :automatic");
				}
			}
			else {
				System.out.println("# of linear iterations:"+nlit);
			}
			System.out.println("# of non-linear iterations:"+nnit);
			System.out.println("concentration map output:"+iconc);
			System.out.println("spherical charge distbn:"+isph);
			System.out.println("INSIGHT format output:"+ibios);
			System.out.println("site potential output:"+isite);
		}
		System.out.println("modified atom file output:"+iatout);
		System.out.println("map file label:");
		System.out.println(toplbl);
		//
		if(ipdbrd||ipdbwrt) {
	  		if(ipdbrd) System.out.println("set to  read unformatted pdb file");
			if(ipdbrd&&ipdbwrt) {
				System.out.println("");
				System.out.println("WARNING: can not write an unformatted pdb");
				System.out.println("file, while reading in an unformatted pdb file");
				System.out.println("Therefore the write option is disabled");
			}else {
				if(ipdbwrt) System.out.println("set to write unformatted pdb file");
			}
		}
		//
		if(ifrcrd||ifrcwrt) {
			if(ifrcrd&&isite) System.out.println("set to read unformatted frc.pdb file");
			if(ifrcrd&&ifrcwrt) {
				System.out.println("");
				System.out.println("WARNING: can not write an unformatted frc");
				System.out.println("file, while reading in an unformatted frc file");
				System.out.println("Therefore the write option is disabled");
			}
			else {
				if(ifrcwrt) System.out.println("set to write unformatted frc.pdb file");
			}
		}
		//
		if(!igraph) System.out.println("convergence graph turned off");
		if(!ipoten) System.out.println("potential listings turned off");
		if((icon1!=10)||(icon2!=1)) {
			System.out.println("convergence test interval is every"+icon1+"loops");
			System.out.println("testing"+100/icon2+"%");
		}
		System.out.println(" ");
		//
	}
 	
	
	void off(){
 
		if(iacent) {
			oldmid[0]=acent[0];
			oldmid[1]=acent[1];
			oldmid[2]=acent[2];
			return;
		}
		 
		for( int i = 1;i<=3;i++) {
	  		oldmid[i-1] = oldmid[i-1] - offset[i-1]/scale;
		}
 		
	}


	void grdatm(){
 
		float rmid=(float)((igrid+1)/2);
 		float cmid[]=oldmid;
		float xn1[]=atpos;
		int j;
		for(int i=1;i<=natom;i++) {
			j=newIndexTwo(i,1,3);
			xn2[j]=(xn1[j] -cmid[0])*scale + rmid;
			xn2[j+1]=(xn1[j+1] -cmid[1])*scale + rmid;
			xn2[j+2]=(xn1[j+2] -cmid[2])*scale + rmid;
			
		}  
	}


	

	void crgarr( ) {
  
		float xn1[]=atpos;
		int ic1;
 
		// atmcrg contians grid positions of all charges AND the charge in the 4th field
 
		int j,k,n,ix;
		ic1=0;
		for( ix=1;ix<=natom;ix++) {
			if(Math.abs(chrgv4[ix-1])>1.e-6) {
				ic1=ic1+1;
				j=newIndexTwo(ic1,1,4);
				k=newIndexTwo(ix,1,3);
				atmcrg[j]=xn2[k];
				atmcrg[j+1]=xn2[k+1];
				atmcrg[j+2]=xn2[k+2];
				n=newIndexTwo(ic1,1,3);
				chgpos[n]=xn1[k];
				chgpos[n+1]=xn1[k+1];
				chgpos[n+2]=xn1[k+2];
				atmcrg[j+3]=chrgv4[ix-1];
				crgatn[ic1-1]=ix;
				//System.out.println("atmcrg"+ix+" "+atmcrg[j]+" "+atmcrg[j+1]+" "+atmcrg[j+2]+" "+chrg);
			}
		}
 
		// assign charges for boundary conditions
 
		// ic1 = number of charges
 
		// find charge moments for dipole approximation
 
		qnet=0.0f;
		qplus=0.0f;
		qmin=0.0f;
		for( ix=1;ix<=3;ix++) {
			cqplus[ix-1]=0.0f;
			cqmin[ix-1]=0.0f;
		}
		float chrg;
		for( ix=1;ix<=ic1;ix++) {
			j=newIndexTwo(ix,1,4);
			chrg=atmcrg[j+3];
			qnet=qnet + chrg;
			if(chrg>0.) {
				qplus=qplus + chrg;
				cqplus[0]=cqplus[0] + chrg*atmcrg[j];
				cqplus[1]=cqplus[1] + chrg*atmcrg[j+1];
				cqplus[2]=cqplus[2] + chrg*atmcrg[j+2];
				//System.out.println("atmcrg"+ix+" "+atmcrg[j]+" "+atmcrg[j+1]+" "+atmcrg[j+2]+" "+chrg);
			}else {
				qmin=qmin + chrg;
				cqmin[0]=cqmin[0] + chrg*atmcrg[j];
				cqmin[1]=cqmin[1] + chrg*atmcrg[j+1];
				cqmin[2]=cqmin[2] + chrg*atmcrg[j+2];
			}
		}
 
		// divide by charge totals
 
		if(qplus>1.e-6) {
	  		for(k = 1;k<=3;k++) {
	    			cqplus[k-1] = cqplus[k-1]/qplus;
	  		}
		}
		if(Math.abs(qmin)>1.e-6) {
	  		for( k = 1;k<=3;k++) {
	    			cqmin[k-1] = cqmin[k-1]/qmin;
          		}
		}
	 
		// select those charges which will be charging the grid
 
	    	nqass=ic1;
	    	int rgrid=igrid;
	    	int ic2=0;

		for( ix=1;ix<=ic1;ix++) {
	    		j=newIndexTwo(ix,1,4);
	    		if((atmcrg[j]>1.)&&(atmcrg[j]<rgrid)) {
	    		if((atmcrg[j+1]>1.)&&(atmcrg[j+1]<rgrid)) {
	    		if((atmcrg[j+2]>1.)&&(atmcrg[j+2]<rgrid)) {
	      			ic2=ic2+1;				
	    		}
	   		}
	    		}
		}
		nqgrd=ic2;
		
		chrgv2=new float[4*nqgrd];
		ic2=0;
		for( ix=1;ix<=ic1;ix++) {
	    		j=newIndexTwo(ix,1,4);
	    		if((atmcrg[j]>1.)&&(atmcrg[j]<rgrid)) {
	    		if((atmcrg[j+1]>1.)&&(atmcrg[j+1]<rgrid)) {
	    		if((atmcrg[j+2]>1.)&&(atmcrg[j+2]<rgrid)) {
	      			ic2=ic2+1;
				k=newIndexTwo(ix,1,4);
	      			chrgv2[k]=atmcrg[j];
	      			chrgv2[k+1]=atmcrg[j+1];
	      			chrgv2[k+2]=atmcrg[j+2];
	      			chrgv2[k+3]=atmcrg[j+3];
	    		}
	   		}
	    		}
		}
		
	}


	
	void wrtadt() {
 
		System.out.println("  ");
		System.out.println("box fill      (%):"+perfil);
		System.out.println("xmin,xmax     (A):"+cmin[0]+" "+cmax[0]);
		System.out.println("ymin,ymax     (A):"+cmin[1]+" "+cmax[1]);
		System.out.println("zmin,zma      (A):"+cmin[2]+" "+cmax[2]);
		System.out.println("x,y,z range   (A):"+cran[0]+" "+cran[1]+" "+cran[2]);
		System.out.println("scale   (grids/A):"+scale);
		System.out.println("object centre (A):"+oldmid[0]+" "+oldmid[1]+" "+oldmid[2]);
		System.out.println("number of atom coordinates read  :"+natom);
		if(isolv){
			System.out.println("total number of charged atoms    :"+nqass);
     			System.out.println("net assigned charge              :"+qnet);
			System.out.println("assigned positive charge         :"+qplus);
			System.out.println("centred at (gu) :                :"+cqplus[0]+" "+cqplus[1]+" "+cqplus[2]);
			System.out.println("assigned negative charge         :"+qmin);
			System.out.println("centred at (gu) :                 "+cqmin[0]+" "+cqmin[1]+" "+cqmin[2]);
		}
		System.out.println("   ");
	}



	
	void epsmak(){

		float xn1[]=atpos;
  		int i,j,k,n;
		bndeps= new int[igrid*igrid*igrid];
 
		if(epsin==epsout) {
			System.out.println("not going to calculate boundary elements since");
			System.out.println("uniform dielectric");
			ibnum=0;
			return;
		}
 
		// lepsx.y.z and uepsx.y.z should be the upper and lower limits of the
		// expanded box. if the molecule is smaller than this { reduce leps 
		// and upeps accordingly
		// note leps/ueps not yet defined..
 
		j=newIndexTwo(1,1,3);
		float xmin=xn2[j];
		float ymin=xn2[j+1];
		float zmin=xn2[j+2];
		float xmax=xn2[j];
		float ymax=xn2[j+1];
		float zmax=xn2[j+2];
		float rmax=rad3[0];
		

		for( i=2;i<=natom;i++) {
			j=newIndexTwo(i,1,3);
			xmin=Math.min(xn2[j],xmin);
			ymin=Math.min(xn2[j+1],ymin);
			zmin=Math.min(xn2[j+2],zmin);
			xmax=Math.max(xn2[j],xmax);
			ymax=Math.max(xn2[j+1],ymax);
			zmax=Math.max(xn2[j+2],zmax);
			rmax=Math.max(rad3[i-1],rmax);
		}	

		if(rionst!=0)rmax=Math.max(rmax,exrad);

		rmax=rmax*scale;

		xmin=xmin-rmax;
		ymin=ymin-rmax;
		zmin=zmin-rmax;
		xmax=xmax+rmax;
		ymax=ymax+rmax;
		zmax=zmax+rmax;
		j=newIndexTwo(1,1,3);
		limeps[j+2]=Math.max((int)(zmin)-2,1);
		limeps[j+1]=Math.max((int)(ymin)-2,1);
		limeps[j]=Math.max((int)(xmin)-2,1);
 
		j=newIndexTwo(2,1,3);
		limeps[j+2]=Math.min((int)(zmax)+2,igrid);
		limeps[j+1]=Math.min((int)(ymax)+2,igrid);
		limeps[j]=Math.min((int)(xmax)+2,igrid);
 
		for( k=1;k<=igrid;k++){
		for( j=1;j<=igrid;j++){
		for( i=1;i<=igrid;i++){
			n=newIndexFour(i,j,k,1,igrid);
			iepsmp[n]=0;
			iepsmp[n+1]=0;
			iepsmp[n+2]=0;
			n=newIndexThree(i,j,k,igrid);
			idebmap[n]=1;
		}
		}
		}
 
		float hgrl=1.f/2.f/scale;
		// if radprb is less than half of grid spacing, { use old algorithm
		// (sri 29March 93)
 
		// The new algorithm should be able to handle all scales; still remains to be
		// tested (Sri Apr 95)
 		long finish0=finish;
	
		finish=cputime(finish);
		 
		String day=getTimeStr();
		System.out.println("start vw surface at: "+day );
		float s=radprb;
		//radprb=0;
		setout(0.0f);
		//System.out.println(radprb);
		//radprb=s;
		finish=cputime(finish);  
		System.out.println("fill in re-entrant regions at:"+finish );
		//vwtms(ibnum,xn1,natom,oldmid,limeps);
		vwtms();
 

		finish=cputime(finish0);
		if(!isrf) System.out.println("time to turn everything in is:"+finish );
 
		// comment out membrane stuff for a moment.. 		 
	}

	void setout(float radprb){
 		float xn[]=new float[3];
		float sq1[]=new float[31];
		float sq2[]=new float[31];
		float sq3[]=new float[31]; 		
		float rad2a1[]=new float[31]; 
		float rad2a2[]=new float[31]; 
		float rad2a3[]=new float[31];
		int ismin[]=new int[3];
		int ismax[]=new int[3];
		boolean itobig,itest2;

 		float temp;
	 	int iv=0,i=0,k=0;		
		float rad,rad5,radp,rad4,rad2,radp2,rad2a;		 		
		int itest1;		
		int num=0;
		int ixn1=0,iyn1=0,izn1=0;
		float fxn1,fxn2,fxn3;
		float rad2ax,rad2ay,rad2az;
		float temp1,temp2,temp3,distsq,distsq3,distsq2,distsq1,dxyz3,dxyz2,dxyz1 ;
		int i1,i2,i3;
		int iboxt=0,ibox=0;
		float radmax2=0.0f;
		int ix=0;
		int igrdc=0,ioff[]=null,iy=0,iz=0,itest=0,je=0;
		float dist,dist1,dist2,dist3,radtest;
		int lim,limmax;
		// a non-zero entry in iepsmp indicates an atom # plus 1 (to properly treat
		// re-entrant mid-points later (15 Aug 93)
		 
		
 
		
		for( ix=1;ix<=natom;ix++) {
			radmax2=Math.max(radmax2,rad3[ix-1]);
		}
		temp=Math.max(exrad,radprb);
		radmax2=scale*(radmax2+temp);
		lim= (int)(1+ radmax2);
 
		limmax = 12;
		itobig=false;
		if(lim>limmax) itobig=true;
 		
		if(!itobig) {
			radtest= (float) Math.pow((radmax2 + 0.5*Math.sqrt(3.0)),2);
			 
			ibox=0;
			igrdc=(int) (Math.pow((2*lim+1),3));
			ioff= new int[3*igrdc];
			
			for( ix=-lim;ix<=lim;ix++) {
	 		for( iy=-lim;iy<=lim;iy++) {
	   	 	for( iz=-lim;iz<=lim;iz++) {
	    			dist=(float) (Math.pow(ix,2) + Math.pow(iy,2) + Math.pow(iz,2));
	    			dist1=dist + ix + 0.25f;
	    			dist2=dist + iy + 0.25f;
	    			dist3=dist + iz + 0.25f;
	    			itest=0;
	    			if(dist<radtest) itest=1;
	    			if(dist1<radtest) itest=1;
	    			if(dist2<radtest) itest=1;
	   			if(dist3<radtest) itest=1;
	    			if(itest==1) {
					//je=newIndexTwo(ibox,1,3);
	    				ibox=ibox+1;
					je=newIndexTwo(ibox,1,3);
	    				ioff[je]=ix;
	    				ioff[je+1]=iy;
	    				ioff[je+2]=iz;
					//System.out.println("ibox:"+ibox+" "+je+" "+ioff[je]+" "+ioff[je+1]+" "+ioff[je+2]);
	    			}
			}
			}
			}
		}
 
		// set interiors
 		float nnn=0;
		for( iv=1; iv<=natom;iv++) {
 
			// restore values
 
			rad= rad3[iv-1];
			i=newIndexTwo(iv,1,3);
			xn[0]=xn2[i];
			xn[1]=xn2[i+1];
			xn[2]=xn2[i+2];
			if(rad<1.e-6) continue;
 
			// scale radius to grid
 
	  		rad = rad*scale;
	  		rad5= (float) (Math.pow((rad + 0.5f),2));
	  		radp = rad + exrad*scale;
	  		rad = rad + radprb*scale;
	  		rad4= (float) (Math.pow((rad + 0.5f),2));
	  		rad2 = rad*rad;
	  		radp2 = radp*radp;
 
			// set dielectric map
 
			// check if sphere sits within limits of box
			itest2=false;
        		for( k = 1;k<=3;k++) {
          			ismin[k-1] = (int) (xn[k-1] - radmax2 - 1.);
	  			itest1=ismin[k-1];
	    			ismin[k-1] = Math.min(ismin[k-1],igrid);
	    			ismin[k-1] = Math.max(ismin[k-1],1);
	    			if(itest1!=ismin[k-1]) itest2=true;
          			ismax[k-1] = (int)(xn[k-1] + radmax2 + 1.);
	  			itest1=ismax[k-1];
	    			ismax[k-1] = Math.min(ismax[k-1],igrid);
	    			ismax[k-1] = Math.max(ismax[k-1],1);
	    			if(itest1!=ismax[k-1]) itest2=true;
			}
 			
			if(itest2||itobig) {
				// slow method
				num=num+1;
				rad2a = rad2 - 0.25f;
          			for( iz =  ismin[2];iz<=ismax[2];iz++){
            			for( iy =  ismin[1];iy<=ismax[1];iy++){
             			for( ix =  ismin[0];ix<=ismax[0];ix++) {
					
                 			dxyz1 = (ix - xn[0]);
                  			dxyz2 = (iy - xn[1]);
                  			dxyz3 = (iz - xn[2]);
                  			distsq =(float)(Math.pow(dxyz1,2) +Math.pow(dxyz2,2) +Math.pow(dxyz3,2));
		  			distsq1 = distsq + dxyz1;
		  			distsq2 = distsq + dxyz2;
		 			distsq3 = distsq + dxyz3;
					i=newIndexFour(ix,iy,iz,1,igrid);
					if(distsq1<rad2a) iepsmp[i]=iv+1;
					if(distsq2<rad2a) iepsmp[i+1]=iv+1;
					if(distsq3<rad2a) iepsmp[i+2]=iv+1;
					i=newIndexThree(ix,iy,iz,igrid);
					if(distsq<radp2)  idebmap[i] =0;
					nnn+=distsq1+distsq2+distsq3;
					//if(nnn<100) {
					//System.out.println(nnn+":"+ix+" "+iy+" "+iz+" "+distsq1+" "+distsq2+" "+distsq3);
					//System.out.println(nnn+":"+ix+" "+iy+" "+iz+" "+iepsmp[i]+" "+iepsmp[i+1]+" "+iepsmp[i+2]);
					//}
				}
				}
				}
			}else {
				// faster method
				rad2a = rad2 - 0.25f;
				
				ixn1=Math.round(xn[0]);
				iyn1=Math.round(xn[1]);
				izn1=Math.round(xn[2]);
				fxn1=ixn1-xn[0];
				fxn2=iyn1-xn[1];
				fxn3=izn1-xn[2];
				rad2ax=rad2a-fxn1;
				rad2ay=rad2a-fxn2;
				rad2az=rad2a-fxn3;
				for(ix=-lim;ix<=lim;ix++) {
					temp1=ix+fxn1;
					temp2=ix+fxn2;
					temp3=ix+fxn3;
					//System.out.println("ix "+ix+" "+temp1+" "+temp2+" "+temp3);
					sq1[ix+15]=temp1*temp1;
					sq2[ix+15]=temp2*temp2;
					sq3[ix+15]=temp3*temp3;
					rad2a1[ix+15]=rad2a-temp1;
					rad2a2[ix+15]=rad2a-temp2;
					rad2a3[ix+15]=rad2a-temp3;
				}
				//$DIR NO_RECURRENCE
		  		for( i=1;i<=ibox;i++) {
		  			je=newIndexTwo(i,1,3);
		  			i1= ioff[je];
		  			i2= ioff[je+1];
		  			i3= ioff[je+2];
		  			ix=ixn1+ i1;
		  			iy=iyn1+ i2;
		  			iz=izn1+ i3;
        				distsq = sq1[i1+15] +sq2[i2+15] + sq3[i3+15];
					je=newIndexFour(ix,iy,iz,1,igrid);
					if(distsq<rad2a1[i1+15]) iepsmp[je]=iv+1;
					if(distsq<rad2a2[i2+15]) iepsmp[je+1]=iv+1;
					if(distsq<rad2a3[i3+15]) iepsmp[je+2]=iv+1;
					je=newIndexThree(ix,iy,iz,igrid);
        				if(distsq<radp2)   idebmap[je]=0;
				}
			}
 
		}
		System.out.println("nnn:"+nnn);
		ioff= null;
 
	}

	void vwtms(){

		
		// subroutine to take a van der Waals epsmap and expand it into a
		// molecular volume map..
		// procedure is to first generate a list of accessible points
		// whose collective volume may contain any midpoint in the box.
		// next send these points to be mapped into indexing arrays...
		// { take the vanderwaals epsmap and grow it out to the molecular
		// epsmap, checking against the set of accessible points..
		//
		// can't do the vw surface by projecting midpoints to accessible surface and
		// testing if that point is in or outside the accessible volume (20 Oct 92)
		// (but can do the contact part of MS that way, Nov 93)
		//
		float xn1[]=atpos;
		float s1,s2,s3,xx,yy,zz,dmn,dist;	
		boolean outcb[]=new boolean[125];  
		String line;
		int nbra[]=new int[1000];
		int iy,je,iz,iv,iall=0,egird[],n1,n2,iarv,iac=0,ia,kk;
		int eps[]=new int[6],ncms=0;;
 		int ix2,iy2,iz2,it1,it2,it3,nnbr=0;
		float x1,x1x,x1y,x1z,r0a,prbrd2,rsm,dis,dr2,dr1,dr3,u1,u2,u3;
		float dx1,dx2,dx3,ds2,dsr;	
		boolean out,intb,exists,goto201;
		boolean nbe[]=new boolean[7];
		
		float xg[]=new float[3];
		float del,cbln,cba;
		int jx,jy,jz,iacv,ii,liml,limu;
		
		int epsn,i,j,k;
		int egrid[];
		float goff[]=new float[18];
		float off;
		int ix,n,nmt,nmmt,nn;
 
		int ibmx=500000;
		r0= new float[natom];
		r02=new float[natom];
		rs2= new float[natom];
		ast= new int[natom];
		//ibnd= new int[3*ibmx];
		
		for( i=-2;i<=2;i++){
		for( j=-2;j<=2;j++){
		for( k=-2;k<=2;k++){
			epsn=newIndexThree(i,j,k,-2,5);
			outcb[epsn]=true;
		}
		}
		}
		for( i=-1;i<=1;i++){
		for( j=-1;j<=1;j++){
		for( k=-1;k<=1;k++){
			epsn=newIndexThree(i,j,k,-2,5);
			outcb[epsn]=false;
		}
		}
		}

		nbe[0]=false;
		nbe[6]=false;
		nbe[1]=true;
		nbe[2]=true;
		nbe[3]=true;
		nbe[4]=true;
		nbe[5]=true;
		//
		
		off=0.5f/scale;
		for(i=1;i<=6;i++){
		for(j=1;j<=3;j++){
			epsn=newIndexTwo(i,j,3);
			goff[epsn]=0.0f;
		}
		}
		epsn=newIndexTwo(1,1,3);
		goff[epsn]=off;
		epsn=newIndexTwo(2,2,3);
		goff[epsn]=off;
		epsn=newIndexTwo(3,3,3); 
		goff[epsn]=off;
		epsn=newIndexTwo(4,1,3); 
		goff[epsn]=-off;
		epsn=newIndexTwo(5,2,3); 
		goff[epsn]=-off;
		epsn=newIndexTwo(6,3,3); 
		goff[epsn]=-off;
	
		// convertion from grid to float coordinates(can also use routine gtoc)
		
		x1=1.0f/scale;
		x1x=oldmid[0]-(1.0f/scale)*(float)(igrid+1)*0.5f;
		x1y=oldmid[1]-(1.0f/scale)*(float)(igrid+1)*0.5f;
		x1z=oldmid[2]-(1.0f/scale)*(float)(igrid+1)*0.5f;
		// find extrema
		cmin[0]=6000f;
		cmin[1]=6000f;
		cmin[2]=6000f;
		cmax[0]=-6000f;
		cmax[1]=-6000f;
		cmax[2]=-6000f;
		rdmx=0.0f;
		//System.out.println("natom:"+natom);
		for(ix=1;ix<=natom;ix++) {
			epsn=newIndexTwo(ix,1,3);
			cmin[0]=Math.min(cmin[0],xn1[epsn]);
			cmin[1]=Math.min(cmin[1],xn1[epsn+1]);
			cmin[2]=Math.min(cmin[2],xn1[epsn+2]);
			cmax[0]=Math.max(cmax[0],xn1[epsn]);
			cmax[1]=Math.max(cmax[1],xn1[epsn+1]);
			cmax[2]=Math.max(cmax[2],xn1[epsn+2]);
			rdmx=Math.max(rdmx,rad3[ix-1]);
			//System.out.println("ix:"+ix+" rdmx:"+rdmx);
		}
		 
		//

		int jje,kke;
		n=0;
		jje=newIndexTwo(1,1,3);
		kke=newIndexTwo(2,1,3);
		for( k=limeps[jje+2]+1;k<=limeps[kke+2]-1;k++)  
		for( j=limeps[jje+1]+1;j<=limeps[kke+1]-1;j++)  
		for( i=limeps[jje]+1;i<=limeps[kke]-1;i++) {	
			nn=0;
			epsn=newIndexFour(i,j,k,1,igrid);
			if(iepsmp[epsn]>0)nn=nn+1;
			if(iepsmp[epsn+1]>0)nn=nn+1;
			if(iepsmp[epsn+2]>0)nn=nn+1;
			epsn=newIndexFour(i-1,j,k,1,igrid);
			if(iepsmp[epsn]>0)nn=nn+1;
			epsn=newIndexFour(i,j-1,k,2,igrid);
			if(iepsmp[epsn]>0)nn=nn+1;
			epsn=newIndexFour(i,j,k-1,3,igrid);
			if(iepsmp[epsn]>0)nn=nn+1;
			//if(nbe[nn]) {
				n=n+1;
			//}
		} 
		System.out.println("ibmx "+n+" "+(3*n));
		ibmx=n+1;
		ibnd= new int[3*ibmx];
		//
		//System.out.println("nnn "+n+" "+(3*n));
		// find vanderwaals boundary
		n=0;
		nmt=0;
		nmmt=0;
		// NB change limits to those of the molecule.
		// set for iepsmp NOT equal to unity
		 
		jje=newIndexTwo(1,1,3);
		kke=newIndexTwo(2,1,3);
		for( k=limeps[jje+2]+1;k<=limeps[kke+2]-1;k++) {
		for( j=limeps[jje+1]+1;j<=limeps[kke+1]-1;j++) {
		for( i=limeps[jje]+1;i<=limeps[kke]-1;i++) {
			nn=0;
			epsn=newIndexFour(i,j,k,1,igrid);
			if(iepsmp[epsn]>0)nn=nn+1;
			if(iepsmp[epsn+1]>0)nn=nn+1;
			if(iepsmp[epsn+2]>0)nn=nn+1;
			epsn=newIndexFour(i-1,j,k,1,igrid);
			if(iepsmp[epsn]>0)nn=nn+1;
			epsn=newIndexFour(i,j-1,k,2,igrid);
			if(iepsmp[epsn]>0)nn=nn+1;
			epsn=newIndexFour(i,j,k-1,3,igrid);
			if(iepsmp[epsn]>0)nn=nn+1;

			if(nbe[nn]) {
				n=n+1;
				epsn=newIndexThree(i,j,k,igrid);
				bndeps[epsn]=n;
				ibnd[3*n-3]=i;
				ibnd[3*n-2]=j;
				ibnd[3*n-1]=k;
			}else {
				epsn=newIndexThree(i,j,k,igrid);
				bndeps[epsn]=0;
			}

		}
		}
		}
		ibnum=n;
		System.out.println("boundary points on vw surface= "+ibnum);
		if(ibnum>ibmx){
			System.out.println("ibnum= "+ibnum+" is greater than ibmx = "+ibmx);
			System.out.println("increase ibmx in vwtms.f");
			return;
		}
		//
		
		for( i=1;i<=natom;i++) {
			r0a=rad3[i-1]+radprb;
			r0[i-1]=r0a;
			r02[i-1]=r0a*r0a;
			rs2[i-1]=0.99999f*r02[i-1];
		}
		// if prbrad=0.0, quit
		if(radprb<1.e-6){
			ibgrd= new int[3*ibnum];
			for( i=1;i<=3*ibnum;i++) {
				ibgrd[i-1]=ibnd[i-1];
			}
			for(i=1;i<=natom;i++) {
				ast[i-1]=0;
			}

			// goto 222;

			if(isolv&&(irea||logs||isen||isch)) {
				System.out.println("scaling boundary grid points");
				scspos=new float[3*ibnum];
				scsnor=new float[3*ibnum];

				for( j=1;j<=ibnum;j++) {
					scspos[3*j-2-1]=(float)(ibgrd[3*j-2-1]);
					scspos[3*j-1-1]=(float)(ibgrd[3*j-1-1]);
					scspos[3*j-1]=(float)(ibgrd[3*j-1]);
				}

				//sclbp(natom,igrid,xn1,scale,radprb,oldmid,ibnum, extot,iall,scspos,scsnor);
				iall=sclbp(xn1,iall,scspos,scsnor);
				//call dtime(tary);
				//System.out.println("time taken = "+tary[1-1]);
				System.out.println(" iall points had to be assigned by global comparison");
			}

			if(isrf){
				//i_egrid = 0;haha
        			egrid= new int[3*igrid*igrid*igrid];
				//**msrf(natom,igrid,xn1,scale,radprb,oldmid,extot,ibnum,egrid);
				msrf(egrid,xn1);
				egrid= null;
        			//i_egrid= memalloc[i_egrid,0];
			}

        		iab1= null;
        		iab2= null;
        		icume= null;

			r0= null;
			r02= null;
			rs2=null;
			ast= null;
 			return;
		}

        	if(iacs) System.out.println("  opening surface file");
	
		// make a list of accessible points..,expos. all scaling of grid
		// points will be done to thses points..
		prbrd2=radprb*radprb;
 

		// calculate an identity for this conformation
		rsm=0.0f;
		for( i=1;i<=natom;i++) {
			epsn=newIndexTwo(i,1,3);
			rsm=rsm+rad3[i-1]*(Math.abs(xn1[epsn])+Math.abs(xn1[epsn+1])+Math.abs(xn1[epsn+2]));
		}
	 
	 
	 	finish=cputime(finish); 
		//System.out.println("rdmx:"+rdmx);
		sas(xn1,natom,radprb,extot);
		finish=cputime(finish); 
		System.out.println("mkacc time = "+finish);

		
		
		del=1.f/scale;
		del=Math.max(del,radprb);
		cbln=rdmx+del;
		cubedata(2.0f,natom,cbln);
        	cbn1= new int[(lcb+1)*(mcb+1)*(ncb+1)];
        	cbn2= new int[(lcb+1)*(mcb+1)*(ncb+1)];
        	cbal= new int[27*natom];
		cube(natom,xn1,rad3,0.0f);
	 
 
		// link the accessible points into iab1 and iab2
 
		//indverdata(radprb,scale);
		System.out.println("vrws,radprb "+radprb);
		indverdata();
		System.out.println("lcb1,mcb1,ncb1,grdi "+lcb1+" "+mcb1+" "+ncb1+" "+grdi);
		cba=1.f/grdi;
        	iab1= new int[(lcb1+1)*(mcb1+1)*(ncb1+1)];
        	iab2= new int[(lcb1+1)*(mcb1+1)*(ncb1+1)];
        	icume=new int[extot];
		
	 	System.out.println("grid for indexing accessible points = "+cba);
 

		//indver(extot);
		indver();
		// write out surface data
		 
        	if(false) {
			for( i=1;i<=extot;i++) {
				je=newIndexTwo(i,1,3);
				xg[1-1]=expos[je];
				xg[2-1]=expos[je+1];
				xg[3-1]=expos[je+2];
				iv=1;
				//System.out.println("surface "+i+" "+expos[je]+" "+expos[je+1]+" "+expos[je+2]);
				//**watput(i,iv,xo,line);
		
			}	 
		}
		// now start the expansion
		// m1= the number of boundary points removed from list
 
		int ncav=0,mpr,ndv,mr,m;
		// START
		n1=1;
		n2=ibnum;
		// m= number of new boundary elements..
		mpr=100000;
		ndv=0;
		System.out.println("ibnum "+ibnum+" "+prbrd2);
		jx=jy=jz=0;dist=0;
		while(true){
		m=0;
		mr=0;
		for( i=n1;i<=n2;i++) {
	
			ix=ibnd[3*i-3];
			iy=ibnd[3*i-2];
			iz=ibnd[3*i-1];
			je=newIndexThree(ix,iy,iz,igrid);
			if(bndeps[je]!=0){ 
				je=newIndexFour(ix,iy,iz,1,igrid);
				eps[0]=iepsmp[je];
				eps[1]=iepsmp[je+1];
				eps[2]=iepsmp[je+2];
				je=newIndexFour(ix-1,iy,iz,1,igrid);
				eps[3]=iepsmp[je];
				je=newIndexFour(ix,iy-1,iz,2,igrid);
				eps[4]=iepsmp[je];
				je=newIndexFour(ix,iy,iz-1,3,igrid);
				eps[5]=iepsmp[je];
 
				xg[0]=(float)(ix)*x1 +x1x;
				xg[1]=(float)(iy)*x1 +x1y;
				xg[2]=(float)(iz)*x1 +x1z;
 				//System.out.println("xg: "+i+" "+(xg[0])+" "+(xg[1])+" "+(xg[2]));
				re200:
				for(j=1;j<=6;j++) { 
					 //System.out.println("dist..."+j+" "+dist+" "+(eps[j-1]));
				    	 if(eps[j-1]==0) {
 
						// add midpoint offset to grid point..
						je=newIndexTwo(j,1,3);
						s1=xg[0]+goff[je];
						s2=xg[1]+goff[je+1];
						s3=xg[2]+goff[je+2];
						// determine if this virgin midpoint is in or out
						xx=(s1-mnx)*grdi;
						yy=(s2-mny)*grdi;
						zz=(s3-mnz)*grdi;
						jx=(int)(xx);
						jy=(int)(yy);
						jz=(int)(zz);

 

						dmn=1000.f;
						iacv=0;
						je=jx*(mcb1+1)*(ncb1+1)+jy*(ncb1+1)+jz;
        					liml=iab1[je];
        					limu=iab2[je];
						//System.out.println("lim..."+j+" "+liml+" "+limu+" "+(eps[j-1])); 
						for(ii=liml;ii<=limu;ii++) {
							//System.out.println("ii "+ii);
							iarv= icume[ii-1];
							je=newIndexTwo(iarv,1,3);
							dist=(float) (Math.pow((s1-expos[je]),2) 
								+Math.pow((s2-expos[je+1]),2) 
								+Math.pow((s3-expos[je+2]),2));	
							//System.out.println("ii..."+j+" "+ii+" "+dist);
							if(dist<prbrd2){
								eps[j-1]=-1;
								//goto 200;
								continue re200; 
							}else if(dist<dmn){
	 							iacv=iarv;
								dmn=dist;
							}
						}
						 
						// -1,0,0
						jx=jx-1;
        					je=jx*(mcb1+1)*(ncb1+1)+jy*(ncb1+1)+jz;
        					liml=iab1[je];
        					limu=iab2[je];
					
						for( ii=liml;ii<=limu;ii++) {
							iarv= icume[ii-1];
							je=newIndexTwo(iarv,1,3);
							dist=(float) (Math.pow((s1-expos[je]),2) 
								+Math.pow((s2-expos[je+1]),2) 
								+Math.pow((s3-expos[je+2]),2));	
							if(dist<prbrd2){
								eps[j-1]=-1;
								//goto 200;
								continue re200; 
							}else if(dist<dmn){
	 							iacv=iarv;
								dmn=dist;
							}
						}
						// 1,0,0
						jx=jx+2;
        					je=jx*(mcb1+1)*(ncb1+1)+jy*(ncb1+1)+jz;
        					liml=iab1[je];
        					limu=iab2[je];
						for( ii=liml;ii<=limu;ii++) {
							iarv= icume[ii-1];
							je=newIndexTwo(iarv,1,3);
							dist=(float) (Math.pow((s1-expos[je]),2) 
								+Math.pow((s2-expos[je+1]),2) 
								+Math.pow((s3-expos[je+2]),2));	
	 
							if(dist<prbrd2){
								eps[j-1]=-1;
								//goto 200;
								continue re200; 
							}else if(dist<dmn){
	 							iacv=iarv;
								dmn=dist;
							}
						}
						// 0,-1,0
						jx=jx-1;
						jy=jy-1;
        					je=jx*(mcb1+1)*(ncb1+1)+jy*(ncb1+1)+jz;
        					liml=iab1[je];
        					limu=iab2[je];

						for( ii=liml;ii<=limu;ii++) {
							iarv= icume[ii-1];
							je=newIndexTwo(iarv,1,3);
							dist=(float) (Math.pow((s1-expos[je]),2) 
								+Math.pow((s2-expos[je+1]),2) 
								+Math.pow((s3-expos[je+2]),2));	
	 
							if(dist<prbrd2){
								eps[j-1]=-1;
								continue re200; 
							}else if(dist<dmn){
	 							iacv=iarv;
								dmn=dist;
							}
 
						}
						// 0,1,0
						jy=jy+2;
       						je=jx*(mcb1+1)*(ncb1+1)+jy*(ncb1+1)+jz;
        					liml=iab1[je];
        					limu=iab2[je];

						for( ii=liml;ii<=limu;ii++) {
							iarv= icume[ii-1];
							je=newIndexTwo(iarv,1,3);
							dist=(float) (Math.pow((s1-expos[je]),2) 
								+Math.pow((s2-expos[je+1]),2) 
								+Math.pow((s3-expos[je+2]),2));	
	 
							if(dist<prbrd2){
								eps[j-1]=-1;
								continue re200; 
							}else if(dist<dmn){
	 							iacv=iarv;
								dmn=dist;
							}
 
						}
						// 0,0,-1
						jy=jy-1;
						jz=jz-1;
        					je=jx*(mcb1+1)*(ncb1+1)+jy*(ncb1+1)+jz;
        					liml=iab1[je];
        					limu=iab2[je];

						for( ii=liml;ii<=limu;ii++) {
							iarv= icume[ii-1];
							je=newIndexTwo(iarv,1,3);
							dist=(float) (Math.pow((s1-expos[je]),2) 
								+Math.pow((s2-expos[je+1]),2) 
								+Math.pow((s3-expos[je+2]),2));	
	 
							if(dist<prbrd2){
								eps[j-1]=-1;
								continue re200; 
							}else if(dist<dmn){
	 							iacv=iarv;
								dmn=dist;
							} 
						}
						// 0,0,1
						jz=jz+2;
        					je=jx*(mcb1+1)*(ncb1+1)+jy*(ncb1+1)+jz;
        					liml=iab1[je];
        					limu=iab2[je];

						for( ii=liml;ii<=limu;ii++) {
							iarv= icume[ii-1];
							je=newIndexTwo(iarv,1,3);
							dist=(float) (Math.pow((s1-expos[je]),2) 
								+Math.pow((s2-expos[je+1]),2) 
								+Math.pow((s3-expos[je+2]),2));	
	 
							if(dist<prbrd2){
								eps[j-1]=-1;
								continue re200; 
							}else if(dist<dmn){
	 							iacv=iarv;
								dmn=dist;
							} 
						}
						// nn=2
						// 1,0,1
						jx=jx+1;
        					je=jx*(mcb1+1)*(ncb1+1)+jy*(ncb1+1)+jz;
        					liml=iab1[je];
        					limu=iab2[je];

	 					for( ii=liml;ii<=limu;ii++) {
		 					iarv= icume[ii-1];
	 						je=newIndexTwo(iarv,1,3);
	 						dist=(float) (Math.pow((s1-expos[je]),2) 
								+Math.pow((s2-expos[je+1]),2) 
								+Math.pow((s3-expos[je+2]),2));	
	 
	 						if(dist<prbrd2){
	 							eps[j-1]=-1;
	 							continue re200; 
	 						}else if(dist<dmn){
	  							iacv=iarv;
	 							dmn=dist;
				 			} 
	 					}
 						// -1,0,1
	 					jx=jx-2;
         					je=jx*(mcb1+1)*(ncb1+1)+jy*(ncb1+1)+jz;
        					liml=iab1[je];
        					limu=iab2[je];

	 					for( ii=liml;ii<=limu;ii++) {
	 						iarv= icume[ii-1];
	 						je=newIndexTwo(iarv,1,3);
	 						dist=(float) (Math.pow((s1-expos[je]),2) 
								+Math.pow((s2-expos[je+1]),2) 
								+Math.pow((s3-expos[je+2]),2));	
	 
	 						if(dist<prbrd2){
	 							eps[j-1]=-1;
	 							continue re200; 
	 						}else if(dist<dmn){
	  							iacv=iarv;
	 							dmn=dist;
	 						} 
	 					}
 						// 0,1,1
	 					jx=jx+1;
	 					jy=jy+1;
         					je=jx*(mcb1+1)*(ncb1+1)+jy*(ncb1+1)+jz;
        					liml=iab1[je];
        					limu=iab2[je];

	 					for( ii=liml;ii<=limu;ii++) {
	 						iarv= icume[ii-1];
	 						je=newIndexTwo(iarv,1,3);
	 						dist=(float) (Math.pow((s1-expos[je]),2) 
								+Math.pow((s2-expos[je+1]),2) 
								+Math.pow((s3-expos[je+2]),2));	
	 
	 						if(dist<prbrd2){
	 							eps[j-1]=-1;
	 							continue re200; 
	 			 			}else if(dist<dmn){
	  							iacv=iarv;
	 							dmn=dist;
	 						} 
	 					}
 						// 0,-1,1
	 					jy=jy-2;
         					je=jx*(mcb1+1)*(ncb1+1)+jy*(ncb1+1)+jz;
        					liml=iab1[je];
        					limu=iab2[je];

	 					for( ii=liml;ii<=limu;ii++) {
	 						iarv= icume[ii-1];
	 						je=newIndexTwo(iarv,1,3);
	 						dist=(float) (Math.pow((s1-expos[je]),2) 
								+Math.pow((s2-expos[je+1]),2) 
								+Math.pow((s3-expos[je+2]),2));	
	 
	 						if(dist<prbrd2){
	 							eps[j-1]=-1;
	 							continue re200; 
	 						}else if(dist<dmn){
	  							iacv=iarv;
	 							dmn=dist;
	 						} 
	 					}
 						// -1,-1,0
	 					jz=jz-1;
	 					jx=jx-1;
         					je=jx*(mcb1+1)*(ncb1+1)+jy*(ncb1+1)+jz;
        					liml=iab1[je];
        					limu=iab2[je];

	 					for( ii=liml;ii<=limu;ii++) {
	 						iarv= icume[ii-1];
	 						je=newIndexTwo(iarv,1,3);
	 						dist=(float) (Math.pow((s1-expos[je]),2) 
								+Math.pow((s2-expos[je+1]),2) 
								+Math.pow((s3-expos[je+2]),2));	
	 
	 						if(dist<prbrd2){
	 							eps[j-1]=-1;
	 							continue re200; 
	 						}else if(dist<dmn){
	  							iacv=iarv;
	 							dmn=dist;
	 						} 
	 					}
 						// 1,-1,0
	 					jx=jx+2;
         					je=jx*(mcb1+1)*(ncb1+1)+jy*(ncb1+1)+jz;
        					liml=iab1[je];
        					limu=iab2[je];

	 					for( ii=liml;ii<=limu;ii++) {
	 						iarv= icume[ii-1];
	 						je=newIndexTwo(iarv,1,3);
	 						dist=(float) (Math.pow((s1-expos[je]),2) 
								+Math.pow((s2-expos[je+1]),2) 
								+Math.pow((s3-expos[je+2]),2));	
	 
	 						if(dist<prbrd2){
	 							eps[j-1]=-1;
	 							continue re200; 
	 						}else if(dist<dmn){
	  							iacv=iarv;
	 							dmn=dist;
	 						} 
	 					}
 						// 1,1,0
	 					jy=jy+2;
         					je=jx*(mcb1+1)*(ncb1+1)+jy*(ncb1+1)+jz;
        					liml=iab1[je];
        					limu=iab2[je];

	 					for( ii=liml;ii<=limu;ii++) {
	 						iarv= icume[ii-1];
	 						je=newIndexTwo(iarv,1,3);
	 						dist=(float) (Math.pow((s1-expos[je]),2) 
								+Math.pow((s2-expos[je+1]),2) 
								+Math.pow((s3-expos[je+2]),2));	
	 
	 						if(dist<prbrd2){
	 							eps[j-1]=-1;
	 							continue re200; 
	 						}else if(dist<dmn){
	  							iacv=iarv;
	 							dmn=dist;
	 						} 
	 					}
 				 		// -1,1,0
	 					jx=jx-2;
         					je=jx*(mcb1+1)*(ncb1+1)+jy*(ncb1+1)+jz;
        					liml=iab1[je];
        					limu=iab2[je];

	 					for( ii=liml;ii<=limu;ii++) {
	 						iarv= icume[ii-1];
	 						je=newIndexTwo(iarv,1,3);
	 						dist=(float) (Math.pow((s1-expos[je]),2) 
								+Math.pow((s2-expos[je+1]),2) 
								+Math.pow((s3-expos[je+2]),2));	
	 
	 						if(dist<prbrd2){
	 							eps[j-1]=-1;
	 							continue re200; 
	 						}else if(dist<dmn){
	  							iacv=iarv;
	 							dmn=dist;
	 						} 
	 					}
 						// -1,0,-1
	 					jz=jz-1;
	 					jy=jy-1;
        					je=jx*(mcb1+1)*(ncb1+1)+jy*(ncb1+1)+jz;
        					liml=iab1[je];
        					limu=iab2[je];

	 					for( ii=liml;ii<=limu;ii++) {
	 						iarv= icume[ii-1];
	 						je=newIndexTwo(iarv,1,3);
	 						dist=(float) (Math.pow((s1-expos[je]),2) 
								+Math.pow((s2-expos[je+1]),2) 
								+Math.pow((s3-expos[je+2]),2));	
	 
	 						if(dist<prbrd2){
	 							eps[j-1]=-1;
	 							continue re200; 
	 						}else if(dist<dmn){
	  							iacv=iarv;
	 				 			dmn=dist;
	 						} 
	 					}
 						// 1,0,-1
	 					jx=jx+2;
         					je=jx*(mcb1+1)*(ncb1+1)+jy*(ncb1+1)+jz;
        					liml=iab1[je];
        					limu=iab2[je];

	 					for( ii=liml;ii<=limu;ii++) {
	 						iarv= icume[ii-1];
	 						je=newIndexTwo(iarv,1,3);
	 						dist=(float) (Math.pow((s1-expos[je]),2) 
								+Math.pow((s2-expos[je+1]),2) 
								+Math.pow((s3-expos[je+2]),2));	
	 
	 						if(dist<prbrd2){
	 							eps[j-1]=-1;
	 							continue re200; 
	 						}else if(dist<dmn){
	  							iacv=iarv;
	 							dmn=dist;
	 						} 
	 					}
 						// 0,1,-1
	 					jx=jx-1;
	 					jy=jy+1;
         					je=jx*(mcb1+1)*(ncb1+1)+jy*(ncb1+1)+jz;
        					liml=iab1[je];
        					limu=iab2[je];

	 					for( ii=liml;ii<=limu;ii++) {
	 						iarv= icume[ii-1];
	 						je=newIndexTwo(iarv,1,3);
	 						dist=(float) (Math.pow((s1-expos[je]),2) 
								+Math.pow((s2-expos[je+1]),2) 
								+Math.pow((s3-expos[je+2]),2));	
	 
	 						if(dist<prbrd2){
	 							eps[j-1]=-1;
	 							continue re200; 
	 						}else if(dist<dmn){
	  							iacv=iarv;
	 							dmn=dist;
	 						} 
	 					}
 						// 0,-1,-1
	 					jy=jy-2;
         					je=jx*(mcb1+1)*(ncb1+1)+jy*(ncb1+1)+jz;
        					liml=iab1[je];
        					limu=iab2[je];

	 					for( ii=liml;ii<=limu;ii++) {
	 						iarv= icume[ii-1];
	 						je=newIndexTwo(iarv,1,3);
	 						dist=(float) (Math.pow((s1-expos[je]),2) 
								+Math.pow((s2-expos[je+1]),2) 
								+Math.pow((s3-expos[je+2]),2));	
	 
	 						if(dist<prbrd2){
	 							eps[j-1]=-1;
	 							continue re200; 
	 						}else if(dist<dmn){
	  							iacv=iarv;
	 							dmn=dist;
	 						} 
	 					}
 						// nn=3
 						// -1,-1,-1
	 					jx=jx-1;
         					je=jx*(mcb1+1)*(ncb1+1)+jy*(ncb1+1)+jz;
        					liml=iab1[je];
        					limu=iab2[je];

	 					for( ii=liml;ii<=limu;ii++) {
	 						iarv= icume[ii-1];
	 						je=newIndexTwo(iarv,1,3);
	 						dist=(float) (Math.pow((s1-expos[je]),2) 
								+Math.pow((s2-expos[je+1]),2) 
								+Math.pow((s3-expos[je+2]),2));	
	 
	 						if(dist<prbrd2){
	 							eps[j-1]=-1;
	 							continue re200; 
	 						}else if(dist<dmn){
	  							iacv=iarv;
	 							dmn=dist;
	 						} 
	 					}
 						// 1,-1,-1
	 					jx=jx+2;
         					je=jx*(mcb1+1)*(ncb1+1)+jy*(ncb1+1)+jz;
        					liml=iab1[je];
        					limu=iab2[je];

	 					for( ii=liml;ii<=limu;ii++) {
		 					iarv= icume[ii-1];
		 					je=newIndexTwo(iarv,1,3);
		 					 
	 						dist=(float) (Math.pow((s1-expos[je]),2) 
								+Math.pow((s2-expos[je+1]),2) 
								+Math.pow((s3-expos[je+2]),2));	
		 					if(dist<prbrd2){
			 					eps[j-1]=-1;
			 					continue re200; 
		 					}else if(dist<dmn){
	 		 					iacv=iarv;
			 					dmn=dist;
		 					} 
	 					}
 						// 1,1,-1
		 				jy=jy+2;
         					je=jx*(mcb1+1)*(ncb1+1)+jy*(ncb1+1)+jz;
        					liml=iab1[je];
        					limu=iab2[je];

	 					for( ii=liml;ii<=limu;ii++) {
	 						iarv= icume[ii-1];
	 						je=newIndexTwo(iarv,1,3);
	 						dist=(float) (Math.pow((s1-expos[je]),2) 
								+Math.pow((s2-expos[je+1]),2) 
								+Math.pow((s3-expos[je+2]),2));		
	 
	 						if(dist<prbrd2){
	 							eps[j-1]=-1;
	 							continue re200; 
	 						}else if(dist<dmn){
	  							iacv=iarv;
	 							dmn=dist;
	 						} 
	 					}
 						// -1,1,-1
	 					jx=jx-2;
         					je=jx*(mcb1+1)*(ncb1+1)+jy*(ncb1+1)+jz;
        					liml=iab1[je];
        					limu=iab2[je];

	 					for( ii=liml;ii<=limu;ii++) {
	 						iarv= icume[ii-1];
	 						je=newIndexTwo(iarv,1,3);
	 						dist=(float) (Math.pow((s1-expos[je]),2) 
								+Math.pow((s2-expos[je+1]),2) 
								+Math.pow((s3-expos[je+2]),2));	
	 
	 						if(dist<prbrd2){
	 							eps[j-1]=-1;
	 							continue re200; 
	 						}else if(dist<dmn){
	  							iacv=iarv;
	 							dmn=dist;
	 						} 
	 					}
 						// -1,1,1
	 					jz=jz+2;
        					je=jx*(mcb1+1)*(ncb1+1)+jy*(ncb1+1)+jz;
        					liml=iab1[je];
        					limu=iab2[je];

	 					for( ii=liml;ii<=limu;ii++) {
	 						iarv= icume[ii-1];
	 						je=newIndexTwo(iarv,1,3);
	 						dist=(float) (Math.pow((s1-expos[je]),2) 
								+Math.pow((s2-expos[je+1]),2) 
								+Math.pow((s3-expos[je+2]),2));	
	 
	 						if(dist<prbrd2){
	 							eps[j-1]=-1;
	 							continue re200; 
	 						}else if(dist<dmn){
	  							iacv=iarv;
	 							dmn=dist;
	 						} 
	 					}
 						// 1,1,1
	 					jx=jx+2;
         					je=jx*(mcb1+1)*(ncb1+1)+jy*(ncb1+1)+jz;
        					liml=iab1[je];
        					limu=iab2[je];

	 					for( ii=liml;ii<=limu;ii++) {
	 						iarv= icume[ii-1];
	 						je=newIndexTwo(iarv,1,3);
	 						dist=(float) (Math.pow((s1-expos[je]),2) 
								+Math.pow((s2-expos[je+1]),2) 
								+Math.pow((s3-expos[je+2]),2));	
	 
	 						if(dist<prbrd2){
	 							eps[j-1]=-1;
	 				 			continue re200; 
	 						}else if(dist<dmn){
	  							iacv=iarv;
	 							dmn=dist;
	 						} 
	 					}
 						// 1,-1,1
	 					jy=jy-2;
         					je=jx*(mcb1+1)*(ncb1+1)+jy*(ncb1+1)+jz;
        					liml=iab1[je];
        					limu=iab2[je];

	 					for( ii=liml;ii<=limu;ii++) {
	 						iarv= icume[ii-1];
	 						je=newIndexTwo(iarv,1,3);
	 						dist=(float) (Math.pow((s1-expos[je]),2) 
								+Math.pow((s2-expos[je+1]),2) 
								+Math.pow((s3-expos[je+2]),2));	
	 
	 						if(dist<prbrd2){
	 							eps[j-1]=-1;
	 							continue re200; 
	 						}else if(dist<dmn){
	  							iacv=iarv;
	 							dmn=dist;
	 						} 
	 					}
 						// -1,-1,1
	 					jx=jx-2;
         					je=jx*(mcb1+1)*(ncb1+1)+jy*(ncb1+1)+jz;
        					liml=iab1[je];
        					limu=iab2[je];
	 					for( ii=liml;ii<=limu;ii++) {
	 						iarv= icume[ii-1];
	 						je=newIndexTwo(iarv,1,3);
	 						dist=(float) (Math.pow((s1-expos[je]),2) 
								+Math.pow((s2-expos[je+1]),2) 
								+Math.pow((s3-expos[je+2]),2));	
	 
	 						if(dist<prbrd2){
	 							eps[j-1]=-1;
	 							continue re200; 
	 						}else if(dist<dmn){
	  							iacv=iarv;
	 							dmn=dist;
	 						} 
	 					}
 
						// it might be in the contact region; find the closest atom surface
						it1=(int)((s1-xo)*cbai);
						it2=(int)((s2-yo)*cbai);
						it3=(int)((s3-zo)*cbai);
						dmn=100.f;
						iac=0;
						nnbr=0;
						 
        					//liml=cbn1[it1+1+(lcb+1)*it2+(lcb+1)*(mcb+1)*it3-1];
        					//limu=cbn2[it1+1+(lcb+1)*it2+(lcb+1)*(mcb+1)*it3-1];
					 	je=it1*(mcb+1)*(ncb+1)+it2*(ncb+1)+it3;
						liml=cbn1[je];
						limu=cbn2[je];
						for( kk=liml;kk<=limu;kk++) {
							ia=cbal[kk-1];
							if(ast[ia-1]==0){
								nnbr=nnbr+1;
								nbra[nnbr-1]=ia;
							}
						} 
						for( ii=1;ii<=nnbr;ii++) {
							ia=nbra[ii-1];
							je=newIndexTwo(ia,1,3);
							dx1=s1-xn1[je];
							dx2=s2-xn1[je+1];
							dx3=s3-xn1[je+2];
							ds2=dx1*dx1+dx2*dx2+dx3*dx3;
							dis=(float)(Math.sqrt(ds2)-rad3[ia-1]);
	 						if(dis<dmn){
	 							dmn=dis;
	 							iac=ia;
	 						}
						}

						if(iac==0){
							//write(6,*)i,j,' might be a cavity point'
	 						ncav=ncav+1;
							//possibly a cavity point
							//goto 201;
						}else {

							// check to see if it is in the contact region of that atom by projecting it to
							// the atom's acc surface and checking against the acc volumes of nearby atoms
							je=newIndexTwo(iac,1,3);
							dr1=s1-xn1[je];
							dr2=s2-xn1[je+1];
							dr3=s3-xn1[je+2];
							dsr=(float) (Math.sqrt(dr1*dr1+dr2*dr2+dr3*dr3));
							u1=xn1[je]+dr1/dsr*r0[iac-1];
							u2=xn1[je+1]+dr2/dsr*r0[iac-1];
							u3=xn1[je+2]+dr3/dsr*r0[iac-1];
							it1=(int)((u1-xo)*cbai);
							it2=(int)((u2-yo)*cbai);
							it3=(int)((u3-zo)*cbai);
							 
        						//liml=cbn1[it1+1+(lcb+1)*it2+(lcb+1)*(mcb+1)*it3-1];
        						//limu=cbn2[it1+1+(lcb+1)*it2+(lcb+1)*(mcb+1)*it3-1];
							je=it1*(mcb+1)*(ncb+1)+it2*(ncb+1)+it3;
							liml=cbn1[je];
							limu=cbn2[je];
							goto201 = false; 
							for( kk=liml;kk<=limu; kk++) {
								ia=cbal[kk-1];
								je=newIndexTwo(ia,1,3);
								dx1=u1-xn1[je];
								dx2=u2-xn1[je+1];
								dx3=u3-xn1[je+2];
								ds2=dx1*dx1+dx2*dx2+dx3*dx3;
								//if(ds2<rs2[ia-1])goto 201;
								if(ds2<rs2[ia-1]) {
									goto201	= true;
									break;
								}
							}
							if(!goto201) {
								// it is in the contact region. flag the midpoint so it is not checked again
								eps[j-1]=-iac;
								continue re200; 
							} 
						}
						
						//201	
						eps[j-1]=1;
						// remap iepsmp
						je=newIndexFour(ix,iy,iz,1,igrid);
        					iepsmp[je]=eps[0];
        					iepsmp[je+1]=eps[1];
        					iepsmp[je+2]=eps[2];
						je=newIndexFour(ix-1,iy,iz,1,igrid);
        					iepsmp[je]=eps[3];
						je=newIndexFour(ix,iy-1,iz,2,igrid);
        					iepsmp[je]=eps[4];
						je=newIndexFour(ix,iy,iz-1,3,igrid);
        					iepsmp[je]=eps[5];
						//
						// check to see if the nearest neighbour status has been changed..
						ix2=ix;
						iy2=iy;
						iz2=iz;
						// if the nearest neighbour is a box boundary point { skip this
						// since box boundary points can not also be dielctric boundary points
						if(j==1) {
							ix2=ix+1;
							if(ix2==igrid)continue re200; 
						}else if(j==2) {
							iy2=iy+1;
							if(iy2==igrid) continue re200; 
						}else if(j==3) {
							iz2=iz+1;
							if(iz2==igrid) continue re200; 
						}else if(j==4) {
							ix2=ix-1;
							if(ix2==1)continue re200; 
						}else if(j==5) {
							iy2=iy-1;
							if(iy2==1) continue re200; 
						}else if(j==6) {
							iz2=iz-1;
							if(iz2==1) continue re200; 
						}
						//
        					nn=0;
						je=newIndexFour(ix2,iy2,iz2,1,igrid);
        					if(iepsmp[je]>0) nn=nn+1;
        					if(iepsmp[je+1]>0) nn=nn+1;
        					if(iepsmp[je+2]>0) nn=nn+1;
						je=newIndexFour(ix2-1,iy2,iz2,1,igrid);
        					if(iepsmp[je]>0) nn=nn+1;
						je=newIndexFour(ix2,iy2-1,iz2,2,igrid);
        					if(iepsmp[je]>0) nn=nn+1;
						je=newIndexFour(ix2,iy2,iz2-1,3,igrid);
        					if(iepsmp[je]>0) nn=nn+1;
						//
						je=newIndexThree(ix2,iy2,iz2,igrid);
						if(nn==6&&bndeps[je]!=0) {
							// reset bndeps for that point.
							bndeps[je]=0;
							mr=mr+1;
							//System.out.println("mr=mr+1 "+j+" "+mr+" "+ix2+" "+iy2+" "+iz2);
						}
						
						if(nn!=6&&bndeps[je]==0) {
							// create a new boundary point..
							m=m+1;
							bndeps[je]=n2+m;
							
							ibnd[3*(n2+m)-3]=ix2;
							ibnd[3*(n2+m)-2]=iy2;
							ibnd[3*(n2+m)-1]=iz2;
						}
						
						//
					}  
				}//200	continue
				//
				// remap iepsmp in case there have been changes..
				// (changes 0's to -1's)
				je=newIndexFour(ix,iy,iz,1,igrid);
        			iepsmp[je]=eps[0];
        			iepsmp[je+1]=eps[1];
        			iepsmp[je+2]=eps[2];
				je=newIndexFour(ix-1,iy,iz,1,igrid);
        			iepsmp[je]=eps[3];
				je=newIndexFour(ix,iy-1,iz,2,igrid);
        			iepsmp[je]=eps[4];
				je=newIndexFour(ix,iy,iz-1,3,igrid);
        			iepsmp[je]=eps[5];
				// check to see if this is still a boundary point
        			nn=0;
				je=newIndexFour(ix,iy,iz,1,igrid);
        			if(iepsmp[je]>0) nn=nn+1;
        			if(iepsmp[je+1]>0) nn=nn+1;
        			if(iepsmp[je+2]>0) nn=nn+1;
				je=newIndexFour(ix-1,iy,iz,1,igrid);
        			if(iepsmp[je]>0) nn=nn+1;
				je=newIndexFour(ix,iy-1,iz,2,igrid);
        			if(iepsmp[je]>0) nn=nn+1;
				je=newIndexFour(ix,iy,iz-1,3,igrid);
        			if(iepsmp[je]>0) nn=nn+1;
				// if not now a boundary element change bndeps
				if(nn==6){
					je=newIndexThree(ix,iy,iz,igrid);
					bndeps[je]=0;
					mr=mr+1;
					//System.out.println("mr=mr+2 "+i+" "+mr+" "+ix+" "+iy+" "+iz);
				}
				//System.out.println("nn "+i+" "+nn);
				//  if end for whether bndeps is nonzero
			}
		//
		// next boundary point
		}
		//
		n1=n2+1;
		n2=n2+m;
		System.out.println("do it ....m= "+m+" mr = "+mr);
		//
       	 	if(m>mpr){
       	 		ndv=ndv+1;
        		if(ndv>2){
        			System.out.println("surface iteration did not converge");
        			return;
        		}
        	}else {
       			 ndv=0;
        	}
		if(m>0) continue;
		else    break;
		}
		if(n2>ibmx){
			System.out.println("ibnd upper bound "+n2+" + exceeds"+ ibmx);
			return;
		}
		//
		cbn1= null;
		cbn2= null;
		cbal=  null;

		finish=cputime(finish);
		System.out.println("time to grow re-entrant surface = "+finish);
		System.out.println("no. cavity mid-points inaccessible to solvent = "+ ncav);

		// consolidate the list, removing dead boundary points, adding new ones..
		j=0;
		ncms=0;
		for( i=1;i<=n2;i++) {
			ix=ibnd[3*i-3];
			iy=ibnd[3*i-2];
			iz=ibnd[3*i-1];
			je=newIndexThree(ix,iy,iz,igrid);
			if(bndeps[je]!=0) {
				j=j+1;				
			}
		}

		if(j>ibmx) {
			System.out.println("no. ms points exceeds ibmx");
			return;
		}
		ibnum=j;
		ibgrd=new int[3*ibnum];

		//ibgrd= new int[3*ibmx];
		j=0;
		ncms=0;
		for( i=1;i<=n2;i++) {
			ix=ibnd[3*i-3];
			iy=ibnd[3*i-2];
			iz=ibnd[3*i-1];
			je=newIndexThree(ix,iy,iz,igrid);
			if(bndeps[je]!=0) {
				j=j+1;
				bndeps[je]=j;
	
				ibgrd[3*j-3]=ix;
				ibgrd[3*j-2]=iy;
				ibgrd[3*j-1]=iz;	 
				
				je=newIndexFour(ix,iy,iz,1,igrid);
				if(iepsmp[je]<0)iepsmp[je]=0;
				if(iepsmp[je+1]<0)iepsmp[je+1]=0;
				if(iepsmp[je+2]<0)iepsmp[je+2]=0;
				je=newIndexFour(ix-1,iy,iz,1,igrid);
				if(iepsmp[je]<0)iepsmp[je]=0;
				je=newIndexFour(ix,iy-1,iz,2,igrid);
				if(iepsmp[je]<0)iepsmp[je]=0;
				je=newIndexFour(ix,iy,iz-1,3,igrid);
				if(iepsmp[je]<0)iepsmp[je]=0;
			}
		}
		
		
		ibnd= null;
		bndeps=null;
		//if(true) return;
		//
		// scale bondary grid point positions relative to acc data
		//222	
		//System.out.println("scaling boundary grid points....");
		if(isolv&&(irea||logs||isen||isch)) {
			System.out.println("scaling boundary grid points");
			scspos=new float[3*ibnum];
			scsnor=new float[3*ibnum];

			for( j=1;j<=ibnum;j++) {
				scspos[3*j-3]=(float)(ibgrd[3*j-3]);
				scspos[3*j-2]=(float)(ibgrd[3*j-2]);
				scspos[3*j-1]=(float)(ibgrd[3*j-1]);
			}

			//sclbp(natom,igrid,xn1,scale,radprb,oldmid,ibnum, extot,iall,scspos,scsnor);
			iall=sclbp(xn1,iall,scspos,scsnor);
			for( j=1;j<=ibnum;j++) {
				//System.out.println("scspos "+j+" "+(scspos[3*j-3])+" "+(scspos[3*j-2])+" "+(scspos[3*j-1]));
 				//System.out.println("scsnor "+i+" "+(scsnor[3*j-3])+" "+(scsnor[3*j-2])+" "+(scsnor[3*j-1]));
			}
			//call dtime(tary);
			finish=cputime(finish);
			System.out.println("time taken = "+finish);
			System.out.println("iall:"+iall);
			System.out.println(" iall points had to be assigned by global comparison");
			//if(true) return;
		}

		if(isrf){
			egrid = null;
        		egrid= new int[3*igrid*igrid*igrid];
			//**msrf(natom,igrid,xn1,scale,radprb,oldmid,extot,ibnum,egrid);
			//msrf(natom,igrid,xn1,scale,radprb,oldmid,extot,ibnum,egrid)
        		//i_egrid= memalloc[i_egrid,0];bbb
			msrf(egrid,xn1);
			egird=null;
		}

        	iab1= null;
        	iab2= null;
        	icume= null;

		r0= null;
		r02= null;
		rs2=null;
		ast= null;
		
	}
	
	
	void indver(){
		 
		// program to compile the lists iab1,iab2 and icume for use in
		// nearest vertex work. iexpos are box coordinates of vertices
		// and dont need to be passed if comaprisons are done with float angstroms..
		// but its often more convenient to do so since grid points have to be
		// converted anyway...
		//
	 
		//int iab1[0:lcb1,0:mcb1,0:ncb1],iab2[0:lcb1,0:mcb1,0:ncb1], icume[1];
		
 		int je,k,j,i,n,ix,iy,iz;
		float x,y,z;
		iexpos= new int[3*extot];
		// initialize grid..
		for( i=0;i<=lcb1;i++) {
		for( j=0;j<=mcb1;j++) {
		for( k=0;k<=ncb1;k++) {
			je=i*(mcb1+1)*(ncb1+1)+j*(ncb1+1)+k;
			iab1[je]=1;
			//je=i*(mcb1+1)*(ncb1+1)+j*(ncb1+1)+k;
			iab2[je]=0;
		}
		}
		}


		
 		System.out.println("lcb1,mcb1,ncb1................. "+extot);
		// make linear arrays for expos
 		System.out.println("lcb1,mcb1,ncb1 "+lcb1+" "+mcb1+" "+ncb1);
		// find the number of points in each box, put in iab2, make iexpos
		for( i=1;i<=extot;i++) {
			
			je=newIndexTwo(i,1,3);
			x=(expos[je]-mnx)*grdi;
			ix=(int)(x);
			y=(expos[je+1]-mny)*grdi;
			iy=(int)(y);
			z=(expos[je+2]-mnz)*grdi;
			iz=(int)(z);
			//System.out.println("je2723:" +i+" "+expos[je]+" "+expos[je+1]+" "+expos[je+2]+" "+mnx+" "+mny+" "+mnz+" "+grdi);
			//System.out.println("expos "+i+": "+(expos[je])+" "+(expos[je+1])+" "+(expos[je+2]));
			je=ix*(mcb1+1)*(ncb1+1)+iy*(ncb1+1)+iz;
			
			//System.out.println("je2723:" +i+" "+extot+" "+ix+" "+iy+" "+iz+" "+je);
			iab2[je]=iab2[je]+1;
			//System.out.println("i,x,y,z "+i+" "+x+" "+y+" "+z);
			je=newIndexTwo(i,1,3);
			iexpos[je]=ix;
			iexpos[je+1]=iy;
			iexpos[je+2]=iz;
		}

		
 
		// check each box for occupancy, using fill number to mark out
		// space in icume, using
		n=0;
		for( i=0;i<=lcb1;i++){
		for( j=0;j<=mcb1;j++){
		for( k=0;k<=ncb1;k++){
 
			// if the box is not empty put start position to n+1 in iab1
			// end to n+box occupancy in iab2, overwriting occupancy..
			je=i*(mcb1+1)*(ncb1+1)+j*(ncb1+1)+k;
			if(iab2[je]!=0) {
				iab1[je]=n+1;
				n=n+iab2[je];
				iab2[je]=n;
				//System.out.println("i,j,k = "+i+" "+j+" "+k+" "+(iab1[je])+" "+(iab2[je]));
			}
 
		}
		}
		}

		
 
		// fill icume using iab1 and iab2, note that iab1 is used to hold the
		// position in icume, therefore needs to be recalculated..

		for( i=1;i<=extot;i++) {
			je=newIndexTwo(i,1,3);
			ix=iexpos[je];
			iy=iexpos[je+1];
			iz=iexpos[je+2];
			je=ix*(mcb1+1)*(ncb1+1)+iy*(ncb1+1)+iz;
			j=iab1[je];
			icume[j-1]=i;
			iab1[je]=iab1[je]+1;
		}
 
		// reset iab1 for use in inner loop
		for( i=1;i<=extot;i++) {
			je=newIndexTwo(i,1,3);
			ix=iexpos[je];
			iy=iexpos[je+1];
			iz=iexpos[je+2];
			je=ix*(mcb1+1)*(ncb1+1)+iy*(ncb1+1)+iz;
			iab1[je]=iab1[je]-1;
		}

 
		for( i=0;i<=lcb1;i++) {
		for( j=0;j<=mcb1;j++) {
		for( k=0;k<=ncb1;k++) {
			je=i*(mcb1+1)*(ncb1+1)+j*(ncb1+1)+k;			
			//System.out.println("i,j,k = "+i+" "+j+" "+k+" "+(iab1[je])+" "+(iab2[je]));
		}
		}
		}
 
		// icume now contains pointers to each dot inside a particular box
		// and each box has 2 pointers into icume.,a start pointer and a finish
		// pointer
		// use however you want...
		//
		iexpos= null;
 
	}

	
	void indverdata(){
	 
		float mxx,mxy,mxz,jig,cba;
		int in,il,idtemp,im; 
		// cba defines scale, cba and jig the box limits..
		// note that cba is the minumum box size. it can be bigger if need be.
		// i.e. can rescale to a larger box if need be.
		float sav=prbrad;
		prbrad=radprb; 
		cba=prbrad;
	
		jig=2.0f*(1.0f/scale);

		mnx=(cmin[0]-rdmx-prbrad);
		mny=(cmin[1]-rdmx-prbrad);
		mnz=(cmin[2]-rdmx-prbrad);
		mxx=cmax[0]+rdmx+prbrad;
		mxy=cmax[1]+rdmx+prbrad;
		mxz=cmax[2]+rdmx+prbrad;
		System.out.println("rdmx,prbrad,cmin[1],cmax[1]"+" "+rdmx+" "+prbrad+" "+(cmin[1])+" "+(cmax[1]));
 		//System.out.println("mxx,mxy,mxz,mnx,mny,mnz"+" "+mxx+" "+mxy+" "+mxz+" "+mnx+
		//		" "+mny+" "+mnz);
		
		// accessible surface points are prbrad away from the actual surface
		// mid points are at most (1/scale) away from the actual surface.
		// midpoints must never be in box zero. accessible point can be in box zero
		// but not in box -1.
		//
		// e.g. mnx+prbrad==van der Waals surface
		// vanderwaals surface - (1/scale)= minumum midpoint position.
		// therefore mnx+prbrad-(1/scale) gt 1
		// and mnx gt 0. the latter is always true since cba is ge. prbrad.
		// therefore we have...
 
		mxx=mxx+jig;
		mxy=mxy+jig;
		mxz=mxz+jig;
		mnx=(mnx-jig);
		mny=(mny-jig);
		mnz=(mnz-jig);
 		//System.out.println("lcb1,mcb1,ncb1,mxx,mxy,mxz,mnx,mny,mnz"+
		//		" "+lcb1+" "+mcb1+" "+ncb1+" "+mxx+" "+mxy+" "+mxz+" "+mnx+
		//		" "+mny+" "+mnz);
		
		re100:
		while(true) {	
			mxx=mxx+cba-prbrad;
			mxy=mxy+cba-prbrad;
			mxz=mxz+cba-prbrad;
			mnx=(mnx-cba+prbrad);
			mny=(mny-cba+prbrad);
			mnz=(mnz-cba+prbrad);
			//
			il=(int)((mxx-mnx)/cba)+1;
			im=(int)((mxy-mny)/cba)+1;
			in=(int)((mxz-mnz)/cba)+1;
 
			// if points are too widely seperated for the scale and idmax {
			// rescale..
			if((il>idmax)||(im>idmax)||(in>idmax)) {
				idtemp=Math.max(il,im);
				idtemp=Math.max(idtemp,in);
				cba=cba*(float)(idtemp+1)/(float)(idmax);
				System.out.println("initial cube size too small, ");
				System.out.println("in assigning accessible points to a grid");
				System.out.println("therefore rescaling...");
				continue re100;
				//**goto 100;
			}
			else break re100;
		}
		lcb1=il;
		mcb1=im;
		ncb1=in;
 		
		// grdi is just the inverse of cba..,used more...
		grdi=1.0f/cba;
		prbrad=sav;
	}



	void mkvtl(int ivbeg,int ivend,int itbeg,int itend,int vtpnt[],int  vtlen[], 
		   int vindx[],int  vtlst[], int tmlst[], int imxtri,int  imxvtx) {
 
		int intsiz;
		int mxtri=300000;
		int vtemp[]=new int[mxtri], i, j1, j2, j3, k1, k2, k3;
		int im, it1, it2, j,je;
		float tarray1[]=new float[2],tarray2[]=new float[2],t1;
		 
 		System.out.println("mkvtl "+ivbeg+" "+ivend+" "+itbeg+" "+itend);
		for( i=ivbeg;i<=ivend;i++) {
	    		vtlen[i-1]=0;
		}

		for( i=itbeg;i<=itend;i++) {
	    		j1=vindx[3*i-3];
	    		j2=vindx[3*i-2];
	    		j3=vindx[3*i-1];
			//System.out.println("j1j2j3 "+i+" "+j1+" "+j2+" "+j3);
	    		vtlen[j1-1]=vtlen[j1-1]+1;
	    		vtlen[j2-1]=vtlen[j2-1]+1;
	    		vtlen[j3-1]=vtlen[j3-1]+1;
		}

		if(ivbeg!=1) vtpnt[ivbeg-1]=vtpnt[ivbeg-1-1]+vtlen[ivbeg-1-1];
		if(ivbeg==1) vtpnt[ivbeg-1]=1;

		for( i=ivbeg+1;i<=ivend;i++){
	    		vtpnt[i-1]=vtpnt[i-1-1]+vtlen[i-1-1];
		}
		vtpnt[ivend+1-1]=vtpnt[ivend-1]+vtlen[ivend-1];

		for( i=ivbeg;i<=ivend;i++) {
	    		vtlen[i-1]=0;
		}
 
		for( i=itbeg;i<=itend;i++) {
	    		j1=vindx[3*i-3];
	    		j2=vindx[3*i-2];
	    		j3=vindx[3*i-1];
	   		k1=vtpnt[j1-1]+vtlen[j1-1];
	    		k2=vtpnt[j2-1]+vtlen[j2-1];
	   		k3=vtpnt[j3-1]+vtlen[j3-1];
	    		vtlen[j1-1]=vtlen[j1-1]+1;
	    		vtlen[j2-1]=vtlen[j2-1]+1;
	    		vtlen[j3-1]=vtlen[j3-1]+1;
	    		vtlst[k1-1]=i;
	    		vtlst[k2-1]=i;
	    		vtlst[k3-1]=i;
		}
	
		//t1=etime[tarray1-1];
 
		for( i=itbeg;i<=itend;i++) {
	    		vtemp[i-1]=0;
		}
 
		im=0;
		for( i=itbeg;i<=itend;i++) {
 
	    		j1=vindx[3*i-3];
	    		j2=vindx[3*i-2];
	    		j3=vindx[3*i-1];
 
			// which triangle borders edge j1-j2
	    		for( j=vtpnt[j1-1];j<=vtpnt[j1-1]+vtlen[j1-1]-1;j++) {
				k1=vtlst[j-1];
				vtemp[k1-1]=j1;
	    		}
	    		it1=0;
	   		for(j=vtpnt[j2-1];j<=vtpnt[j2-1]+vtlen[j2-1]-1;j++) {
				k1=vtlst[j-1];
				if((vtemp[k1-1]==j1)&&(k1!=i)) it1=k1;
				vtemp[k1-1]=j2;
	    		}
			// which point completes the triangle
	    		it2=0;
	    		if(it1!=0) {
				k1=vindx[3*it1-2-1];
				k2=vindx[3*it1-1-1];
				k3=vindx[3*it1-1];
				if((k1!=j1)&&(k1!=j2)) it2=k1;
				if((k2!=j1)&&(k2!=j2)) it2=k2;
				if((k3!=j1)&&(k3!=j2)) it2=k3;
	    		}
			// fill tmlst
	    		je=newIndexTwo(i,1,9);
	    		tmlst[je]=it1;
	    		tmlst[je+3]=it2;
	    		tmlst[je+6]=j3;
	    		im=im+it1;
 
			// which triangle borders edge j2-j3
	    		it1=0;
	    		for( j=vtpnt[j3-1];j<=vtpnt[j3-1]+vtlen[j3-1]-1;j++) {
				k1=vtlst[j-1];
				if((vtemp[k1-1]==j2)&&(k1!=i)) it1=k1;
				vtemp[k1-1]=j3;
	    		}
			// which point completes the triangle
	    		it2=0;
	    		if(it1!=0) {
				k1=vindx[3*it1-3];
				k2=vindx[3*it1-2];
				k3=vindx[3*it1-1];
				if((k1!=j2)&&(k1!=j3)) it2=k1;
				if((k2!=j2)&&(k2!=j3)) it2=k2;
				if((k3!=j2)&&(k3!=j3)) it2=k3;
	    		}
			// fill tmlst
		 	je=newIndexTwo(i,2,9);
	    		tmlst[je]=it1;
	    		tmlst[je+3]=it2;
	    		tmlst[je+6]=j1;
	    		im=im+it1;
 
			// which triangle borders edge j3-j1
	    		it1=0;
	    		for( j=vtpnt[j1-1];j<=vtpnt[j1-1]+vtlen[j1-1]-1;j++) {
				k1=vtlst[j-1];
				if((vtemp[k1-1]==j3)&&(k1!=i)) it1=k1;
	    		}
			// which point completes the triangle
	    		it2=0;
	    		if(it1!=0) {
				k1=vindx[3*it1-3];
				k2=vindx[3*it1-2];
				k3=vindx[3*it1-1];
				if((k1!=j3)&&(k1!=j1)) it2=k1;
				if((k2!=j3)&&(k2!=j1)) it2=k2;
				if((k3!=j3)&&(k3!=j1)) it2=k3;
	    		}
			// fill tmlst
			je=newIndexTwo(i,3,9);
	    		tmlst[je]=it1;
	    		tmlst[je+3]=it2;
	    		tmlst[je+6]=j2;
	    		im=im+it1;
			//
		}
		 
		//t2=etime[tarray2-1];
 
	}

	
	
	void msrf(int egrid[],float xn1[]) {
		 
		//float xn1[]=atpos;
		float areas,areac,arear,area,zn,yn;
		float xn,rad2,cc,tne4,ss,bb,aa,rad;	 
		//int itot=0,vtot=0;		 
		int je,ib,i,k,j,ia1,ia2,ia3,iate; 
		int imxvtx,mxtri,mxvtx,ntot2=0,imxtri,je2,je1,je3,iv3,iv2,iv1,it;		 
		String fixyn;
		float xo2[]=new float[3],vnorm2[],vmg,vz,vy,vx,tar;
		int isrfin[],iall=0,ntot=0;
		boolean fix;
 		boolean out,nbe[]=new boolean [7],exists; 
		float v1[]=new float[3],v2[]=new float[3],v3[]=new float[3];

		//end of hole-fixing variables
 
		mxvtx=(int)(ibnum*2.2);
		mxtri=2*mxvtx;
		for( i=1;i<=igrid;i++) {
		for( j=1;j<=igrid;j++) {
		for( k=1;k<=igrid;k++) {
			je=newIndexFour(i,j,k,1,igrid);
			egrid[je]=iepsmp[je];
			egrid[je+1]=iepsmp[je+1];
			egrid[je+2]=iepsmp[je+2];
		}
		}
		}
 
		for( i = 1; i<=igrid;i++) {
	 	for( j = 1; j<=igrid;j++) {
			je=newIndexFour(i,1,j,1,igrid);
	   		egrid[je] = egrid[je+1];
			je=newIndexFour(i,j,1,1,igrid);
	    		egrid[je] = egrid[je+2];
	   	}
		}
 
	 	for( k = 2; k<=igrid;k++) {
	 	for( j = 2; j<=igrid;j++) {
		for(i = 2; i<=igrid;i++) {
	    		iate = 0;
			je=newIndexFour(i,j,k,1,igrid);
	    		if (egrid[je] > 0) iate = iate + 1;
	    		if (egrid[je+1]> 0) iate = iate + 1;
	    		if (egrid[je+2] > 0) iate = iate + 1;
			je=newIndexFour(i-1,j,k,1,igrid);
	    		if (egrid[je] > 0) iate = iate + 1;
			je=newIndexFour(i,j-1,k,2,igrid);
	    		if (egrid[je] > 0) iate = iate + 1;
			je=newIndexFour(i,j,k-1,3,igrid);
	    		if (egrid[je] > 0) iate = iate + 1;
	    		if (iate  <= 3) {
				je=newIndexFour(i,j,k,1,igrid);
				egrid[je] = 0;
	    		}else {
				je=newIndexFour(i,j,k,1,igrid);
				egrid[je] = 1;
	    		}
	 	}
	 	}
		}
 		int iii=0;
		for( k = 1; k<=igrid;k++) {
	 	for( j = 1; j<=igrid;j++) {
		for(i = 1; i<=igrid;i++) {
			je=newIndexFour(i,j,k,1,igrid);
			iii+=egrid[je]+egrid[je+1]*9+egrid[je+2];
		}}}

		System.out.println("egrid:"+iii);
			

		//123 	continue

		vindx = new int[mxtri*3];
		vert = new float[mxvtx*3];

	 	//**ex(igrid, egrid, vtot, itot, vindx, vert, "./", 2);
		ex(egrid,vtot,itot,vindx,vert,"./",2);
		  
		//if(true) return;
		if(vtot>mxvtx){
			System.out.println("vtot = "+vtot+" > mxvtx = "+mxvtx);
			System.out.println("increase mxvtx in msrf.f");
			return;
		}

		for(ib=1;ib<=vtot;ib++) {
			je=newIndexTwo(ib,1,3);
			vert[je]=vert[je]/2.f;
			vert[je+1]=vert[je+1]/2.f;
			vert[je+2]=vert[je+2]/2.f;
		}
		//
		itot = itot/3;
		// scale boundary grid point positions relative to acc data
		System.out.println("scaling vertices");
		vnorm = new float[3*vtot];
	 	vnorm2 = null;
		vnorm2 =new float[3*vtot];
		//call dtime(tary);
		 
        	//**sclbp(natom,igrid,xn1,scale,radprb,oldmid,vtot,extot,iall,vert,vnorm);
		for(i=1;i<=vtot;i++) {
			je=newIndexTwo(i,1,3);
			//System.out.println("vert1.."+i+" "+vert[je]+" "+vert[je+1]+" "+vert[je+2]);
			//System.out.println("vnorm1.."+i+" "+vnorm[je]+" "+vnorm[je+1]+" "+vnorm[je+2]);
			//System.out.println("vindx..."+i+" "+vindx[je]+" "+vindx[je+1]+" "+vindx[je+2])
		}
		int itb=ibnum;
		ibnum=vtot;
		iall=sclbp(xn1,iall,vert,vnorm);
		vtot=ibnum;
		ibnum=itb;
		for(i=1;i<=vtot;i++) {
			je=newIndexTwo(i,1,3);
			//System.out.println("myvert.."+i+" "+vert[je]+" "+vert[je+1]+" "+vert[je+2]);
			//System.out.println("vnorm.."+i+" "+vnorm[je]+" "+vnorm[je+1]+" "+vnorm[je+2]);
		}
		
		// fix holes and make vertex to triangle arrays
		// allocate hole-fixing arrays next
		// hole-fixing variables

		if (vtot < mxvtx/2) {
	    		imxvtx = (int)(vtot*2);
		}else {
	    		imxvtx = mxvtx;
		}

		if (itot < mxtri/2) {
	    		imxtri = itot*2;
		}else {
	    		imxtri = mxtri;
		}
	
		vtlen =new int[imxvtx];
		vtlst = new int[6*imxvtx];
		tmlst = new int[9*imxtri];
		vtpnt = new int[imxvtx];
		
         	//**mkvtl(1,vtot, 1,itot, vtpnt, vtlen, vindx, vtlst, tmlst, imxtri, imxvtx);
		mkvtl(1,vtot, 1,itot, vtpnt, vtlen, vindx, vtlst, tmlst, imxtri, imxvtx);
		
		//if(true) return;
        	 
		ntot=fxhl(1,vtot,1,itot,ntot, vtpnt, vtlen, vindx, itot, vtlst, vert, vtot, tmlst, imxvtx, imxtri, vnorm);
		
        	ntot2=fxhl(1,vtot,1,itot,ntot2, vtpnt, vtlen, vindx, itot, vtlst, vert, vtot, tmlst, imxvtx, imxtri, vnorm);
		for(i=1;i<=vtot;i++) {
			je=newIndexTwo(i,1,3);
			//System.out.println("vindx.."+i+" "+vindx[je]+" "+vindx[je+1]+" "+vindx[je+2]+" "+itot+" "+ntot2);
		}
		for(i=1;i<=vtot;i++) {
			je=newIndexTwo(i,1,3);
			//System.out.println("vert.."+i+" "+vert[je]+" "+vert[je+1]+" "+vert[je+2]+" "+itot+" "+ntot2);
		}
 		//if(true) return;

		if (ntot2 > 0)  {
	    		fix = true;
		}else {
	    		fix = false;
		}
		
		//**do while ( fix) {
		while ( fix) {
            		ntot2=fxhl(1,vtot,1,itot,ntot2, vtpnt, vtlen, vindx, itot, vtlst, vert, vtot, tmlst, imxvtx, imxtri, vnorm);
 
	    		if (ntot2 > 0) {
				fix = true;
	    		}else {
				fix = false;
	   		}
		}
		if(itot>mxtri){
			System.out.println("itot = "+itot+" > mxtri = "+mxtri);
			System.out.println("increase mxtri in msrf.f");
			return;
		}
		vtemp = null;

		vtlen=null;
		vtlst=null;
		tmlst=null;
		vtpnt=null;

		System.out.println("number of vertices = "+vtot);
		System.out.println("number of triangles = "+itot);

		for( i=1;i<=vtot;i++) {
			je=newIndexTwo(i,1,3);
			vnorm2[je]=0.0f;
			vnorm2[je+1]=0.0f;
			vnorm2[je+2]=0.0f;
		}
		// calculate area
		area=0.0f;
		areas=0.0f;
		areac=0.0f;
		arear=0.0f;
	
		for( it=1;it<=itot;it++) {
			iv1=vindx[3*it-3];
			iv2=vindx[3*it-2];
			iv3=vindx[3*it-1];
 
			for( k=1;k<=3;k++){
	 			
				v1[k-1]=vert[(iv2-1)*3+k-1]-vert[(iv1-1)*3+k-1];
				v2[k-1]=vert[(iv3-1)*3+k-1]-vert[(iv1-1)*3+k-1];
			}
			vx=v1[2-1]*v2[3-1]-v1[3-1]*v2[2-1];
			vy=v1[3-1]*v2[1-1]-v1[1-1]*v2[3-1];
			vz=v1[1-1]*v2[2-1]-v1[2-1]*v2[1-1];
			vmg=(float)Math.sqrt(vx*vx+vy*vy+vz*vz);
			tar=vmg/2.f;
			je1=newIndexTwo(iv1,1,3);
			je2=newIndexTwo(iv2,1,3);
			je3=newIndexTwo(iv3,1,3);
			vx=vnorm[je1]+vnorm[je2]+vnorm[je3];
			vy=vnorm[je1+1]+vnorm[je2+1]+vnorm[je3+1];
			vz=vnorm[je1+2]+vnorm[je2+2]+vnorm[je3+2];
			vmg=(float) Math.sqrt(vx*vx+vy*vy+vz*vz);
	
			vnorm2[je1]=vnorm2[je1]+vx/vmg;
			vnorm2[je1+1]=vnorm2[je1+1]+vy/vmg;
			vnorm2[je1+2]=vnorm2[je1+2]+vz/vmg;
			vnorm2[je2]=vnorm2[je2]+vx/vmg;
			vnorm2[je2+1]=vnorm2[je2+1]+vy/vmg;
			vnorm2[je2+2]=vnorm2[je2+2]+vz/vmg;
			vnorm2[je3]=vnorm2[je3]+vx/vmg;
			vnorm2[je3+1]=vnorm2[je3+1]+vy/vmg;
			vnorm2[je3+2]=vnorm2[je3+2]+vz/vmg;
			// calculate spherical triangle area if appropriate
			ia1=atndx[iv1-1];
			ia2=atndx[iv2-1];
			ia3=atndx[iv3-1];
			if(ia1>0){
				if(ia1==ia2&&ia1==ia3){
					rad=rad3[ia1-1];
					rad2=rad*rad;
					aa=0.0f;
					bb=0.0f;
					cc=0.0f;
					for( k=1;k<=3;k++) {
						je1=newIndexTwo(iv1,k,3);
						je2=newIndexTwo(iv2,k,3);
						je3=newIndexTwo(iv3,k,3);
						aa=aa+(float)Math.pow((vert[je2]-vert[je1]),2);
						bb=bb+(float)Math.pow((vert[je3]-vert[je2]),2);
						cc=cc+(float)Math.pow((vert[je1]-vert[je3]),2);
					}
					aa=(float)Math.acos(1.-aa/(2.*rad2));
					bb=(float)Math.acos(1.-bb/(2.*rad2));
					cc=(float)Math.acos(1.-cc/(2.*rad2));
					ss=(aa+bb+cc)*.5f;
					//System.out.println("itss "+it+" "+aa+" "+bb+" "+cc);
					tne4=(float)Math.sqrt(Math.tan(ss*.5)*Math.tan((ss-aa)*.5)
					           *Math.tan((ss-bb)*.5)*Math.tan((ss-cc)*.5));
					tar=4.f*(float)Math.atan(tne4)*rad2;
					//System.out.println("itss "+it+" "+aa+" "+bb+" "+cc+" "+tar);
				}
			}
			area=area+tar;
			//System.out.println("area "+it+" "+tar+" "+area);
			//System.out.println("vx "+it+" "+vx+" "+vy+" "+vz);
		}
	
		for( i=1;i<=vtot;i++) {
			je=newIndexTwo(i,1,3);
			xn=vnorm2[je];
			yn=vnorm2[je+1];
			zn=vnorm2[je+2];
			vmg=(float)Math.sqrt(xn*xn+yn*yn+zn*zn);
			je=newIndexTwo(i,1,3);
			vnorm2[je]=vnorm2[je]/vmg;
			vnorm2[je+1]=vnorm2[je+1]/vmg;
			vnorm2[je+2]=vnorm2[je+2]/vmg;
		}

		System.out.println("MS area                = "+area);

	 	//**wrtsurf("grasp.srf",9, oldmid, vert, vtot, vindx, itot, vnorm);
		if(surface==null) surface=new Surface();
		surface.wrtsurf(vtot,itot,vert,vindx,vnorm);
		surface.molecule=molecule;
		surface.radprb=radprb;
		surface.assignAtom();
		//call wrtspdb("test.surf",9,vtot, vert)
	}


	void  cubedata(float fac,int natm,float cbln){
		
		float off=0.1f;
		xo=cmin[0]-fac*cbln-off;
		yo=cmin[1]-fac*cbln-off;
		zo=cmin[2]-fac*cbln-off;
		float xp=cmax[0]+fac*cbln+off;
		float yp=cmax[1]+fac*cbln+off;
		float zp=cmax[2]+fac*cbln+off;
		float bl1=xp-xo;
		float bl2=yp-yo;
		float bl3=zp-zo;
		lcb=(int) (bl1/cbln);
		mcb=(int) (bl2/cbln);
		ncb=(int) (bl3/cbln);
		cbai=1.f/cbln;
	}
	int sclbp(float xn1[],int iall,float vcrd[],float vnr[]) {
		 
		//float xn1[]=atpos;
		//float vcrd[]=scspos;
         	//float vnr[]=scsnor;
		float dr1,dr2,cbln,dsr,dr3,u1,u2,u3,cba,s1,s2,s3;
		float xx,yy,zz,dx,dy,dz,dmn1,dmn2,dmn3,dmx,dcr,ctf,rmn,dist,rdist ;
		int nbra[]=new int[1000],je,ke,i,nnbr,ix,iy,iz,j,ncbp,k;
		int jx,jy,jz,iacl,jjx,jjy,jjz,jxi,jyi,jzi;		  
		boolean out,outcb[]=new boolean[125];
		float xg1,xg2,xg3,dmn,dx1,dx2,dx3,ds2,dis;
		int it1,it2,it3,iac,liml,limu,kk,ia,ii;
		 
		atsurf= new int[ibnum];
		atndx=  new int[ibnum];
		for(i=-2;i<=2;i++) {
		for(j=-2;j<=2;j++) {
		for(k=-2;k<=2;k++) {
			je=newIndexThree(i,j,k,-2,5);
			outcb[je]=true;
		}
		}
		}
		for(i=-1;i<=1;i++) {
		for(j=-1;j<=1;j++) {
		for(k=-1;k<=1;k++) {
			je=newIndexThree(i,j,k,-2,5);
			outcb[je]=false;
		}
		}
		}
		// convertion from grid to float coordinates(can also use routine gtoc)
		float x1=1.0f/scale;
		float x1x=oldmid[0]-(1.0f/scale)*(float)(igrid+1)*0.5f;
		float x1y=oldmid[1]-(1.0f/scale)*(float)(igrid+1)*0.5f;
		float x1z=oldmid[2]-(1.0f/scale)*(float)(igrid+1)*0.5f;
 		
		if(extot==0&&radprb>0.0){
			// find extrema
			cmin[0]=6000;
			cmin[1]=6000;
			cmin[2]=6000;
			cmax[0]=-6000;
			cmax[1]=-6000;
			cmax[2]=-6000;
			rdmx=0.0f;
			for( ix=1;ix<=natom;ix++) {
				je=newIndexTwo(ix,1,3);
				cmin[0]=Math.min(cmin[0],xn1[je]);
				cmin[1]=Math.min(cmin[1],xn1[je+1]);
				cmin[2]=Math.min(cmin[2],xn1[je+2]);
				cmax[0]=Math.max(cmax[0],xn1[je]);
				cmax[1]=Math.max(cmax[1],xn1[je+1]);
				cmax[2]=Math.max(cmax[2],xn1[je+2]);
				rdmx=Math.max(rdmx,rad3[ix-1]);
			}
			//sas(xn1,natom,radprb,extot);
	 		//sas(xn1,natom,radprb,extot,expos,ast);
			sas(xn1,natom,radprb,extot);
		}

		float del=radprb;
		del=Math.max(del,1.f/(2.f*scale));
        	cbln=rdmx+del;
		cubedata(2.0f,natom,cbln);
		System.out.println("lcb.."+lcb+" "+mcb+" "+ncb+" "+rdmx+" "+del);
        	cbn1= new int [(lcb+1)*(mcb+1)*(ncb+1)];
        	cbn2= new int [(lcb+1)*(mcb+1)*(ncb+1)];
        	cbal= new int [27*natom];
		
		//**cube(natom,xn1,rad3,0.0);
		cube(natom,xn1,rad3,0.0f);
	
		ncbp=0;
		
		for(i=1;i<=ibnum;i++) {
			je=newIndexTwo(i,1,3);
        		xg1=vcrd[je]*x1+x1x;
        		xg2=vcrd[je+1]*x1+x1y;
       			xg3=vcrd[je+2]*x1+x1z;
			// find the closest surface atom to the gridpoint
			it1=(int)((xg1-xo)*cbai);
			it2=(int)((xg2-yo)*cbai);
			it3=(int)((xg3-zo)*cbai);
			dmn=100.f;
			iac=0;
			nnbr=0;
 
        		//liml=cbn1[it1+1+(lcb+1)*it2+(lcb+1)*(mcb+1)*it3];
        		//limu=cbn2[it1+1+(lcb+1)*it2+(lcb+1)*(mcb+1)*it3];
			je=it1*(mcb+1)*(ncb+1)+it2*(ncb+1)+it3;
			liml=cbn1[je];
			limu=cbn2[je];
 			//System.out.println("it1,it2,it3 "+i+" "+it1+" "+it2+" "+it3);
			//System.out.println("xg1,xg2,xg3 "+i+" "+xg1+" "+xg2+" "+xg3+" "+liml+" "+limu);
			for( kk=liml;kk<=limu;kk++) {
				ia=cbal[kk-1];
				if(ast[ia-1]==0){
					nnbr=nnbr+1;
					nbra[nnbr-1]=ia;
				}
			}
			for(ii=1;ii<=nnbr;ii++) {
				ia=nbra[ii-1];
				je=newIndexTwo(ia,1,3);
				dx1=xg1-xn1[je];
				dx2=xg2-xn1[je+1];
				dx3=xg3-xn1[je+2];
				ds2=dx1*dx1+dx2*dx2+dx3*dx3;
				dis=(float) (Math.sqrt(ds2)-rad3[ia-1]);
				if(dis<dmn){
					dmn=dis;
					iac=ia;
					
				}
			}
			atsurf[i-1]=iac;
			
	 		if(iac==0){
				System.out.println("no close atom for boundary point "+i);
	 			return iall;
	 		}	
			je=newIndexTwo(iac,1,3);
			dr1=xg1-xn1[je];
			dr2=xg2-xn1[je+1];
			dr3=xg3-xn1[je+2];
			dsr=(float)(Math.sqrt(dr1*dr1+dr2*dr2+dr3*dr3));
			out=true;
			if(radprb>0.0){
				je=newIndexTwo(iac,1,3);
				u1=xn1[je]+dr1/dsr*r0[iac-1];
				u2=xn1[je+1]+dr2/dsr*r0[iac-1];
				u3=xn1[je+2]+dr3/dsr*r0[iac-1];
				it1=(int)((u1-xo)*cbai);
				it2=(int)((u2-yo)*cbai);
				it3=(int)((u3-zo)*cbai);

				nnbr=0;
 
        			//liml=cbn1[it1+1+(lcb+1)*it2+(lcb+1)*(mcb+1)*it3-1];
        			//limu=cbn2[it1+1+(lcb+1)*it2+(lcb+1)*(mcb+1)*it3-1];
				je=it1*(mcb+1)*(ncb+1)+it2*(ncb+1)+it3;
				liml=cbn1[je];
				limu=cbn2[je];
				for(kk=liml;kk<=limu;kk++) {
					ia=cbal[kk-1];
 
					je=newIndexTwo(ia,1,3);
					dx1=u1-xn1[je];
					dx2=u2-xn1[je+1];
					dx3=u3-xn1[je+2];
					ds2=dx1*dx1+dx2*dx2+dx3*dx3;
					if(ds2<rs2[ia-1])out=false;
				}
			}

			if(out){
				ncbp=ncbp+1;
				je=newIndexTwo(i,1,3);
				ke=newIndexTwo(iac,1,3);
        			vcrd[je]=xn1[ke]+dr1*rad3[iac-1]/dsr;
        			vcrd[je+1]=xn1[ke+1]+dr2*rad3[iac-1]/dsr;
       				vcrd[je+2]=xn1[ke+2]+dr3*rad3[iac-1]/dsr;
				je=newIndexTwo(i,1,3);
       				vnr[je]=dr1/dsr;
        			vnr[je+1]=dr2/dsr;
        			vnr[je+2]=dr3/dsr;
				atndx[i-1]=iac;
			}else {
				atndx[i-1]=0;
			}
		}
        	cbn1= null;
        	cbn2= null;
        	cbal= null;

		// scale the re-entrant points with respect to expos
		// if radprb = 0.0 we are done.
		if(radprb>0.0){
			iall=0;

			cba=1.f/grdi;
			for( i=1;i<=ibnum;i++) {
				if(atndx[i-1]==0){
					je=newIndexTwo(i,1,3);
        				s1=vcrd[je]*x1+x1x;
        				s2=vcrd[je+1]*x1+x1y;
        				s3=vcrd[je+2]*x1+x1z;

					xx=(s1-mnx)*grdi;
					yy=(s2-mny)*grdi;
					zz=(s3-mnz)*grdi;
					jx=(int)(xx);
					jy=(int)(yy);
					jz=(int)(zz);
					dx=xx-(float)(jx);
					dy=yy-(float)(jy);
					dz=zz-(float)(jz);
	
					dmn1=Math.min(dx,dy);
					dmn1=Math.min(dmn1,dz);
					dmx=Math.max(dx,dy);
					dmx=Math.max(dmx,dz);
					dmn2=1.0f-dmx;
					dcr=Math.min(dmn1,dmn2);
					ctf=cba*(1+dcr);
					ctf=ctf*ctf;
					iacl=0;
					rmn=100.f;
					
					for( jjx=jx-1;jjx<=jx+1;jjx++){ 
					for( jjy=jy-1;jjy<=jy+1;jjy++){
					for( jjz=jz-1;jjz<=jz+1;jjz++){
					int jje=jjx*(mcb1+1)*(ncb1+1)+jjy*(ncb1+1)+jjz;
					for( ii=iab1[jje];ii<=iab2[jje];ii++) {
						iac= icume[ii-1];
						je=newIndexTwo(iac,1,3);
						dist=(float)(Math.pow((s1-expos[je]),2) +
							Math.pow((s2-expos[je+1]),2) +
							Math.pow((s3-expos[je+2]),2));
						if(dist<rmn){
							rmn=dist;
							iacl=iac;
						}
					}
					}
					}
					}
					if(!(iacl>0&&rmn<ctf))  {
						for(jxi=-2;jxi<=2;jxi++) {
						for(jyi=-2;jyi<=2;jyi++) {
						for(jzi=-2;jzi<=2;jzi++) {
							je=newIndexThree(jxi,jyi,jzi,-2,5);
							if(outcb[je]){
							jjx=jx+jxi;
							if(jjx>=0&&jjx<=lcb1){
							jjy=jy+jyi;
							if(jjy>=0&&jjy<=mcb1){
							jjz=jz+jzi;
							if(jjz>=0&&jjz<=ncb1){
	 							int jje1=jjx*(mcb1+1)*(ncb1+1)+jjy*(ncb1+1)+jjz;
								for(ii=iab1[jje1];ii<=iab2[jje1];ii++) {
									iac= icume[ii-1];
									je=newIndexTwo(iac,1,3);  
									dist=(float) (Math.pow((s1-expos[je]),2) +
									Math.pow((s2-expos[je+1]),2)+
									Math.pow((s3-expos[je+2]),2));
									if(dist<rmn){
										rmn=dist;
										iacl=iac;
									}
								}
							}
							}
							}
							}
						}
						}
						}
						if(!(iacl>0)){
							iall=iall+1;

							for(iac=1;iac<=extot;iac++) {
								je=newIndexTwo(iac,1,3);  
								dist=(float) (Math.pow((s1-expos[je]),2) +
									Math.pow((s2-expos[je+1]),2) +
									Math.pow((s3-expos[je+2]),2));
								if(dist<rmn){
									rmn=dist;
									iacl=iac;
								}
							}
						}
					}
					//300	continue
					//System.out.println("iacl "+iacl);
					je=newIndexTwo(iacl,1,3); 
					
        				dx=s1-expos[je];
        				dy=s2-expos[je+1];
        				dz=s3-expos[je+2];
					rdist=(float) (Math.sqrt(dx*dx+dy*dy+dz*dz));
					if(rdist==0) {
						dist=0.0f;
					}else {
        					dist=radprb/rdist;
					}
					je=newIndexTwo(iacl,1,3); 
					ke=newIndexTwo(i,1,3); 
       			 		vcrd[ke]=expos[je]+dx*dist;
        		 		vcrd[ke+1]=expos[je+1]+dy*dist;
        		 		vcrd[ke+2]=expos[je+2]+dz*dist;
			 		
        		 		if(rdist>1.0e-8){
						ke=newIndexTwo(i,1,3); 
        		 			vnr[ke]=-dx/rdist;
        		 			vnr[ke+1]=-dy/rdist;
        		 			vnr[ke+2]=-dz/rdist;
        		 		}else {
        		 			System.out.println("bdp close to arcp "+i+" "+rdist);
        		 		}

			 	}
			}
		}
		
		System.out.println("% of contact boundary points = "+ ((float)(ncbp)/(float)(ibnum)*100.f));
		
		return iall;
	}


	void cube(int natm,float crd[],float rda[],float radp){
		//include 'acc2.h';
		//float crd[3,natm],rda[natm];
		//int cbn1[0:lcb,0:mcb,0:ncb],cbn2[0:lcb,0:mcb,0:ncb],cbal[1];
		float x,y,z;
		int i,j,k,jx,jy,jz,icum,iz,iy,ix,je;
		icbn=new int[3*natm];

		for( i=0;i<=lcb;i++){
		for( j=0;j<=mcb;j++){
		for( k=0;k<=ncb;k++){
			je=i*(mcb+1)*(ncb+1)+j*(ncb+1)+k;
			cbn1[je]=1;
			cbn2[je]=0;
		}
		}
		}

		for( i=1;i<=natm;i++){
        		if(rda[i-1]>radp){
				je=newIndexTwo(i,1,3);
				x=(crd[je]-xo)*cbai;
				ix=(int)(x);
				y=(crd[je+1]-yo)*cbai;
				iy=(int)(y);
				z=(crd[je+2]-zo)*cbai;
				iz=(int)(z);
	 			if(ix<1||iy<1||iz<1)System.out.println("ix,iy,iz: "+ix+" "+iy+" "+iz);
	 			if(ix>=lcb||iy>=mcb||iz>=ncb)System.out.println("ix,iy,iz: "+ix+" "+iy+" "+iz);
				for( jz=iz-1;jz<=iz+1;jz++) {
	 			for( jy=iy-1;jy<=iy+1;jy++) {
	  			for( jx=ix-1;jx<=ix+1;jx++) {
	  				je=jx*(mcb+1)*(ncb+1)+jy*(ncb+1)+jz;
	   				cbn2[je]=cbn2[je]+1;
					
	  			}
	 			}
				}
				je=newIndexTwo(i,1,3);
        			icbn[je]=ix;
        			icbn[je+1]=iy;
        			icbn[je+2]=iz;
			}
			
		}

		icum=1;
		for( iz=0;iz<=ncb;iz++){
		for( iy=0;iy<=mcb;iy++){
		for( ix=0;ix<=lcb;ix++){
	    		je=ix*(mcb+1)*(ncb+1)+iy*(ncb+1)+iz;
	   		if(cbn2[je]>0){
	   			cbn1[je]=icum;
	   			icum=icum+cbn2[je];
	   		}
	  	}
		}
		}
		if(cbal==null) {
			System.out.println("cbn1 == null");
		}
		for( i=1;i<=natm;i++) {
       		 	if(rda[i-1]>radp){
				je=newIndexTwo(i,1,3);
				ix=icbn[je];
				iy=icbn[je+1];
				iz=icbn[je+2];
	  			je=ix*(mcb+1)*(ncb+1)+iy*(ncb+1)+iz;
				cbal[cbn1[je]-1]=i;
				cbn1[je]=cbn1[je]+1;
			}
		}
		// -1,0,0
		for( i=1;i<=natm;i++){
        		if(rda[i-1]>radp){
				je=newIndexTwo(i,1,3);
				ix=icbn[je];
				iy=icbn[je+1];
				iz=icbn[je+2];
				ix=ix-1;
				je=ix*(mcb+1)*(ncb+1)+iy*(ncb+1)+iz;
				cbal[cbn1[je]-1]=i;
				cbn1[je]=cbn1[je]+1;
				je=newIndexTwo(i,1,3);
				icbn[je]=ix;
				icbn[je+1]=iy;
				icbn[je+2]=iz;
			}
		}
		// 1,0,0
		for( i=1;i<=natm;i++) {
        		if(rda[i-1]>radp){
				je=newIndexTwo(i,1,3);
				ix=icbn[je];
				iy=icbn[je+1];
				iz=icbn[je+2];
				ix=ix+2;
				je=ix*(mcb+1)*(ncb+1)+iy*(ncb+1)+iz;
				cbal[cbn1[je]-1]=i;
				cbn1[je]=cbn1[je]+1;
				je=newIndexTwo(i,1,3);
				icbn[je]=ix;
				icbn[je+1]=iy;
				icbn[je+2]=iz;
			}
		}
		// 0,-1,0
		for( i=1;i<=natm;i++) {
       			if(rda[i-1]>radp){
				je=newIndexTwo(i,1,3);
				ix=icbn[je];
				iy=icbn[je+1];
				iz=icbn[je+2];
				ix=ix-1;
				iy=iy-1;
				je=ix*(mcb+1)*(ncb+1)+iy*(ncb+1)+iz;
				cbal[cbn1[je]-1]=i;
	
				cbn1[je]=cbn1[je]+1;
				je=newIndexTwo(i,1,3);
				icbn[je]=ix;
				icbn[je+1]=iy;
				icbn[je+2]=iz;
			}
		}
		// 0,1,0
		for( i=1;i<=natm;i++) {
        		if(rda[i-1]>radp){
				je=newIndexTwo(i,1,3);
				ix=icbn[je];
				iy=icbn[je+1];
				iz=icbn[je+2];
	
				iy=iy+2;
				je=ix*(mcb+1)*(ncb+1)+iy*(ncb+1)+iz;
				cbal[cbn1[je]-1]=i;
				cbn1[je]=cbn1[je]+1;
				je=newIndexTwo(i,1,3);
				icbn[je]=ix;
				icbn[je+1]=iy;
				icbn[je+2]=iz;
			}
		}
		// 0,0,-1
		for( i=1;i<=natm;i++) {
        		if(rda[i-1]>radp){
				je=newIndexTwo(i,1,3);
				ix=icbn[je];
				iy=icbn[je+1];
				iz=icbn[je+2];
				iy=iy-1;
				iz=iz-1;
				je=ix*(mcb+1)*(ncb+1)+iy*(ncb+1)+iz;
				cbal[cbn1[je]-1]=i;
				cbn1[je]=cbn1[je]+1;
				je=newIndexTwo(i,1,3);
				icbn[je]=ix;
				icbn[je+1]=iy;
				icbn[je+2]=iz;
			}
		}
		// 0,0,1
		for( i=1;i<=natm;i++) {
        		if(rda[i-1]>radp){
				je=newIndexTwo(i,1,3);
				ix=icbn[je];
				iy=icbn[je+1];
				iz=icbn[je+2];
				iz=iz+2;

				je=ix*(mcb+1)*(ncb+1)+iy*(ncb+1)+iz;
				cbal[cbn1[je]-1]=i;
				cbn1[je]=cbn1[je]+1;
				je=newIndexTwo(i,1,3);
				icbn[je]=ix;
				icbn[je+1]=iy;
				icbn[je+2]=iz;	

	 
			}
		}
		// nn=2
		// 1,0,1
		for( i=1;i<=natm;i++) {
        		if(rda[i-1]>radp){
				je=newIndexTwo(i,1,3);
				ix=icbn[je];
				iy=icbn[je+1];
				iz=icbn[je+2];
				ix=ix+1;
	
				je=ix*(mcb+1)*(ncb+1)+iy*(ncb+1)+iz;
				cbal[cbn1[je]-1]=i;
				cbn1[je]=cbn1[je]+1;
				je=newIndexTwo(i,1,3);
				icbn[je]=ix;
				icbn[je+1]=iy;
				icbn[je+2]=iz;	
			}
		}
		// -1,0,1
		for( i=1;i<=natm;i++) {
        		if(rda[i-1]>radp){
				je=newIndexTwo(i,1,3);
				ix=icbn[je];
				iy=icbn[je+1];
				iz=icbn[je+2];
				ix=ix-2;
				je=ix*(mcb+1)*(ncb+1)+iy*(ncb+1)+iz;
				cbal[cbn1[je]-1]=i;
				cbn1[je]=cbn1[je]+1;
				je=newIndexTwo(i,1,3);
				icbn[je]=ix;
				icbn[je+1]=iy;
				icbn[je+2]=iz;	
			}
		}
		// 0,1,1
		for( i=1;i<=natm;i++) {
        		if(rda[i-1]>radp){
				je=newIndexTwo(i,1,3);
				ix=icbn[je];
				iy=icbn[je+1];
				iz=icbn[je+2];
				ix=ix+1;
				iy=iy+1;
				je=ix*(mcb+1)*(ncb+1)+iy*(ncb+1)+iz;
				cbal[cbn1[je]-1]=i;
				cbn1[je]=cbn1[je]+1;
				je=newIndexTwo(i,1,3);
				icbn[je]=ix;
				icbn[je+1]=iy;
				icbn[je+2]=iz;	
			}
		}
		// 0,-1,1
		for( i=1;i<=natm;i++) {
        		if(rda[i-1]>radp){
				je=newIndexTwo(i,1,3);
				ix=icbn[je];
				iy=icbn[je+1];
				iz=icbn[je+2];
				iy=iy-2;
				je=ix*(mcb+1)*(ncb+1)+iy*(ncb+1)+iz;
				cbal[cbn1[je]-1]=i;
				cbn1[je]=cbn1[je]+1;
				je=newIndexTwo(i,1,3);
				icbn[je]=ix;
				icbn[je+1]=iy;
				icbn[je+2]=iz;	
			}
		}
		// -1,-1,0
		for( i=1;i<=natm;i++) {
        		if(rda[i-1]>radp){
				je=newIndexTwo(i,1,3);
				ix=icbn[je];
				iy=icbn[je+1];
				iz=icbn[je+2];
				ix=ix-1;
				iz=iz-1;
				je=ix*(mcb+1)*(ncb+1)+iy*(ncb+1)+iz;
				//System.out.println(ix+" "+iy+" "+iz);
				//System.out.println(i+" "+je);//+" "+(cbn1[je]-1));
				cbal[cbn1[je]-1]=i;
				
				cbn1[je]=cbn1[je]+1;
				je=newIndexTwo(i,1,3);
				icbn[je]=ix;
				icbn[je+1]=iy;
				icbn[je+2]=iz;	
			}
		}
		// 1,-1,0
		for( i=1;i<=natm;i++) {
        		if(rda[i-1]>radp){
				je=newIndexTwo(i,1,3);
				ix=icbn[je];
				iy=icbn[je+1];
				iz=icbn[je+2];
				ix=ix+2;
				je=ix*(mcb+1)*(ncb+1)+iy*(ncb+1)+iz;
				cbal[cbn1[je]-1]=i;
				cbn1[je]=cbn1[je]+1;
				je=newIndexTwo(i,1,3);
				icbn[je]=ix;
				icbn[je+1]=iy;
				icbn[je+2]=iz;	
			}
		}
		// 1,1,0
		for( i=1;i<=natm;i++) {
        		if(rda[i-1]>radp){
				je=newIndexTwo(i,1,3);
				ix=icbn[je];
				iy=icbn[je+1];
				iz=icbn[je+2];
				iy=iy+2;
				je=ix*(mcb+1)*(ncb+1)+iy*(ncb+1)+iz;
				cbal[cbn1[je]-1]=i;
				cbn1[je]=cbn1[je]+1;
				je=newIndexTwo(i,1,3);
				icbn[je]=ix;
				icbn[je+1]=iy;
				icbn[je+2]=iz;	
			}
		}
		// -1,1,0
		for( i=1;i<=natm;i++) {
        		if(rda[i-1]>radp){
				je=newIndexTwo(i,1,3);
				ix=icbn[je];
				iy=icbn[je+1];
				iz=icbn[je+2];
				ix=ix-2;
				je=ix*(mcb+1)*(ncb+1)+iy*(ncb+1)+iz;
				cbal[cbn1[je]-1]=i;
				cbn1[je]=cbn1[je]+1;
				je=newIndexTwo(i,1,3);
				icbn[je]=ix;
				icbn[je+1]=iy;
				icbn[je+2]=iz;	
			}
		}
		// -1,0,-1
		for( i=1;i<=natm;i++) {
        		if(rda[i-1]>radp){
				je=newIndexTwo(i,1,3);
				ix=icbn[je];
				iy=icbn[je+1];
				iz=icbn[je+2];
				iz=iz-1;
				iy=iy-1;
				je=ix*(mcb+1)*(ncb+1)+iy*(ncb+1)+iz;
				cbal[cbn1[je]-1]=i;
				cbn1[je]=cbn1[je]+1;
				je=newIndexTwo(i,1,3);
				icbn[je]=ix;
				icbn[je+1]=iy;
				icbn[je+2]=iz;	
			}
		}
		// 1,0,-1
		for( i=1;i<=natm;i++) {
        		if(rda[i-1]>radp){
				je=newIndexTwo(i,1,3);
				ix=icbn[je];
				iy=icbn[je+1];
				iz=icbn[je+2];
				ix=ix+2;
				je=ix*(mcb+1)*(ncb+1)+iy*(ncb+1)+iz;
				cbal[cbn1[je]-1]=i;
				cbn1[je]=cbn1[je]+1;
				je=newIndexTwo(i,1,3);
				icbn[je]=ix;
				icbn[je+1]=iy;
				icbn[je+2]=iz;	
			}
		}
		// 0,1,-1
		for( i=1;i<=natm;i++) {
        		if(rda[i-1]>radp){
				je=newIndexTwo(i,1,3);
				ix=icbn[je];
				iy=icbn[je+1];
				iz=icbn[je+2];
				ix=ix-1;
				iy=iy+1;
				je=ix*(mcb+1)*(ncb+1)+iy*(ncb+1)+iz;
				cbal[cbn1[je]-1]=i;
				cbn1[je]=cbn1[je]+1;
				je=newIndexTwo(i,1,3);
				icbn[je]=ix;
				icbn[je+1]=iy;
				icbn[je+2]=iz;	
			}
		}
		// 0,-1,-1
		for( i=1;i<=natm;i++) {
        		if(rda[i-1]>radp){
				je=newIndexTwo(i,1,3);
				ix=icbn[je];
				iy=icbn[je+1];
				iz=icbn[je+2];
				iy=iy-2;
				je=ix*(mcb+1)*(ncb+1)+iy*(ncb+1)+iz;
				cbal[cbn1[je]-1]=i;
				cbn1[je]=cbn1[je]+1;
				je=newIndexTwo(i,1,3);
				icbn[je]=ix;
				icbn[je+1]=iy;
				icbn[je+2]=iz;	
			}
		}
		// nn=3
		// -1,-1,-1
		for( i=1;i<=natm;i++) {
        		if(rda[i-1]>radp){
				je=newIndexTwo(i,1,3);
				ix=icbn[je];
				iy=icbn[je+1];
				iz=icbn[je+2];
				ix=ix-1;
				je=ix*(mcb+1)*(ncb+1)+iy*(ncb+1)+iz;
				cbal[cbn1[je]-1]=i;
				cbn1[je]=cbn1[je]+1;
				je=newIndexTwo(i,1,3);
				icbn[je]=ix;
				icbn[je+1]=iy;
				icbn[je+2]=iz;	
			}
		}
		// 1,-1,-1
		for( i=1;i<=natm;i++) {
        		if(rda[i-1]>radp){
				je=newIndexTwo(i,1,3);
				ix=icbn[je];
				iy=icbn[je+1];
				iz=icbn[je+2];
				ix=ix+2;
				je=ix*(mcb+1)*(ncb+1)+iy*(ncb+1)+iz;
				cbal[cbn1[je]-1]=i;
				cbn1[je]=cbn1[je]+1;
				je=newIndexTwo(i,1,3);
				icbn[je]=ix;
				icbn[je+1]=iy;
				icbn[je+2]=iz;	
			}
		}
		// 1,1,-1
		for( i=1;i<=natm;i++) {
        		if(rda[i-1]>radp){
				je=newIndexTwo(i,1,3);
				ix=icbn[je];
				iy=icbn[je+1];
				iz=icbn[je+2];
				iy=iy+2;
				je=ix*(mcb+1)*(ncb+1)+iy*(ncb+1)+iz;
				cbal[cbn1[je]-1]=i;
				cbn1[je]=cbn1[je]+1;
				je=newIndexTwo(i,1,3);
				icbn[je]=ix;
				icbn[je+1]=iy;
				icbn[je+2]=iz;	
			}
		}
		// -1,1,-1
		for( i=1;i<=natm;i++) {
        		if(rda[i-1]>radp){
				je=newIndexTwo(i,1,3);
				ix=icbn[je];
				iy=icbn[je+1];
				iz=icbn[je+2];
				ix=ix-2;
				je=ix*(mcb+1)*(ncb+1)+iy*(ncb+1)+iz;
				cbal[cbn1[je]-1]=i;
				cbn1[je]=cbn1[je]+1;
				je=newIndexTwo(i,1,3);
				icbn[je]=ix;
				icbn[je+1]=iy;
				icbn[je+2]=iz;	
			}
		}
		// -1,1,1
		for( i=1;i<=natm;i++) {
        		if(rda[i-1]>radp){
				je=newIndexTwo(i,1,3);
				ix=icbn[je];
				iy=icbn[je+1];
				iz=icbn[je+2];
				iz=iz+2;
				je=ix*(mcb+1)*(ncb+1)+iy*(ncb+1)+iz;
				cbal[cbn1[je]-1]=i;
				cbn1[je]=cbn1[je]+1;
				je=newIndexTwo(i,1,3);
				icbn[je]=ix;
				icbn[je+1]=iy;
				icbn[je+2]=iz;	
			}
		}
		// 1,1,1
		for( i=1;i<=natm;i++) {
        		if(rda[i-1]>radp){
				je=newIndexTwo(i,1,3);
				ix=icbn[je];
				iy=icbn[je+1];
				iz=icbn[je+2];
				ix=ix+2;
				je=ix*(mcb+1)*(ncb+1)+iy*(ncb+1)+iz;
				cbal[cbn1[je]-1]=i;
				cbn1[je]=cbn1[je]+1;
				je=newIndexTwo(i,1,3);
				icbn[je]=ix;
				icbn[je+1]=iy;
				icbn[je+2]=iz;	
			}
		}
		// 1,-1,1
		for( i=1;i<=natm;i++) {
        		if(rda[i-1]>radp){
				je=newIndexTwo(i,1,3);
				ix=icbn[je];
				iy=icbn[je+1];
				iz=icbn[je+2];
				iy=iy-2;
				je=ix*(mcb+1)*(ncb+1)+iy*(ncb+1)+iz;
				cbal[cbn1[je]-1]=i;
				cbn1[je]=cbn1[je]+1;
				je=newIndexTwo(i,1,3);
				icbn[je]=ix;
				icbn[je+1]=iy;
				icbn[je+2]=iz;	
			}
		}
		// -1,-1,1
		for( i=1;i<=natm;i++) {
        		if(rda[i-1]>radp){
				je=newIndexTwo(i,1,3);
				ix=icbn[je];
				iy=icbn[je+1];
				iz=icbn[je+2];
				ix=ix-2;
				je=ix*(mcb+1)*(ncb+1)+iy*(ncb+1)+iz;
				cbal[cbn1[je]-1]=i;
				cbn1[je]=cbn1[je]+1;
				je=newIndexTwo(i,1,3);
				icbn[je]=ix;
				icbn[je+1]=iy;
				icbn[je+2]=iz;	
			}
		}

		// reset cbn1
		icum=1;
		System.out.println("ncb,mcb,lcb"+ncb+" "+mcb+" "+lcb);
		for(iz=0;iz<=ncb;iz++){
		for(iy=0;iy<=mcb;iy++){
	  	for(ix=0;ix<=lcb;ix++){
	   		je=ix*(mcb+1)*(ncb+1)+iy*(ncb+1)+iz;
			//System.out.println("icum:"+ix+" "+iy+" "+iz+" "+icum+" "+(cbn2[je]));
	   		if(cbn2[je]>0){
	   			cbn1[je]=icum;
	   			icum=icum+cbn2[je];
	   			cbn2[je]=icum-1;
				
	   		}
	  	}
	 	}
		}
		icum=icum-1;
		icbn= null;
		for( i=0;i<=lcb;i++) 
		for( j=0;j<=mcb;j++) 
		for( k=0;k<=ncb;k++) 
		{
			je=i*(mcb+1)*(ncb+1)+j*(ncb+1)+k;
			ix=cbn1[je];
			iy=cbn2[je];
			//System.out.println("cbn1,cbn2 "+i+" "+j+" "+k+" "+ix+" "+iy);
		}
     }



	
	void sas(float crd[],int natm,float radprb,int nacc){
		 
		int mxpr=250000;
		float pi = 3.14159265f;
		int nver=520,nedge=1040;
        	 		 
		int edgv[]=new int[2*nedge],edg[]=new int[nedge];
		int oti[]=new int[nver],st[]=new int[nedge];
		int ic3,ic1,ic2,iv,nbv,nvo,nprx,ii,ip,nxa;
		int je,iv1,iv2,ie,ie2,ie1,ne,ilvl,nlvl,nvi,nv,i,j,nprt,nacct;
		int npr,jj,limu,liml,ix3,ix2,ix1,ke,nst,ia2,k,ia1,nprp;
		float cbln,dctf2,del,d2,dx1,dx2,dx3,dctf,radj,ctf,ctf2,rad,dy2,sm2;
		float xm,rdn,tta,vmg,ym,x3,x2,x1,tm,rv2,cst,rv1,snt,tij2,tij1,sm1;
		float csp,dmg,rvmg,rij,pre,tij3,cf3,cf2,cf1,dy3,ds2,dy1;
		float ver[]=new float[3*nver],rm[]=new float[9];
		nacct=0;
		nprt=0;
		nlvl=3;
		nvi=12;
		pls=new int[2];
		int pls00[];
		cbln=2.f*(rdmx+radprb);
		//System.out.println(natm+" "+cbln+" "+radprb+" "+rdmx);
		cubedata(1.0f,natm,cbln);
		//System.out.println(lcb+" "+mcb+" "+ncb);
      		cbn1= new int[(lcb+1)*(mcb+1)*(ncb+1)];
        	cbn2= new int[(lcb+1)*(mcb+1)*(ncb+1)];
        	cbal= new int[27*natm];
		
		cube(natm,crd,r0,radprb);
		
		
 
		tta=2.f*pi/(float)(nvi);
		for( i=1;i<=nvi;i++) {
			rdn=(i-1)*tta;
			je=newIndexTwo(i,1,3);
			ver[je]=(float)Math.cos(rdn);
			ver[je+1]=(float)Math.sin(rdn);
			ver[je+2]=0.0f;
			j=i+1;
			if(i==nvi)j=1;
			je=newIndexTwo(i,1,2);
			edgv[je]=i;
			edgv[je+1]=j;
		}
		nv=nvi;
		ne=nvi;
		ie1=1;
		ie2=0;

		for( ilvl=1;ilvl<=nlvl;ilvl++) {
			ie1=ie2+1;
			ie2=ne;
			for( ie=ie1;ie<=ie2;ie++) {
				je=newIndexTwo(ie,1,2);
				iv1=edgv[je];
				iv2=edgv[je+1];
				je=newIndexTwo(iv1,1,3);
				xm=ver[je]+ver[newIndexTwo(iv2,1,3)];
				ym=ver[je+1]+ver[newIndexTwo(iv2,2,3)];
				vmg=(float)Math.sqrt(xm*xm+ym*ym);
				nv=nv+1;
				je=newIndexTwo(nv,1,3);
				ver[je]=xm/vmg;
				ver[je+1]=ym/vmg;
				ver[je+2]=0.0f;
				ne=ne+1;
				edg[ie-1]=ne;
				je=newIndexTwo(ne,1,2);
				edgv[je]=iv1;
				edgv[je+1]=nv;
				ne=ne+1;
				je=newIndexTwo(ne,1,2);
				edgv[je]=nv;
				edgv[je+1]=iv2;
			}
		}
		ne=ie2;
		for( ie=ie1;ie<=ne;ie++) {
			edg[ie-1]=-1;
		}
		System.out.println("nv = "+nv+" ne = "+ne);

		for( i=1;i<=natm;i++) {
 
			ast[i-1]=1;
		}

		nacc=0;
		npr=0;
		nprp=0;
		System.out.println("begin....cube");
		int plsnn=0;
		for( i=1;i<=natm;i++) {
			rad=r0[i-1];
			if(rad==radprb) continue;
			je=newIndexTwo(i,1,3);
			x1=crd[je];
			x2=crd[je+1];
			x3=crd[je+2];
			ix1=(int)((x1-xo)*cbai);
			ix2=(int)((x2-yo)*cbai);
			ix3=(int)((x3-zo)*cbai);
 			je=ix1*(mcb+1)*(ncb+1)+ix2*(ncb+1)+ix3;
			//liml=cbn1[ix1+1+(lcb+1)*ix2+(lcb+1)*(mcb+1)*ix3-1];
			//limu=cbn2[ix1+1+(lcb+1)*ix2+(lcb+1)*(mcb+1)*ix3-1];
			//System.out.println("1: liml,limu="+liml+" "+limu);
			//liml=cbn1[ix1+1+(lcb+1)*ix2+(lcb+1)*(mcb+1)*ix3];
			//limu=cbn2[ix1+1+(lcb+1)*ix2+(lcb+1)*(mcb+1)*ix3];
			//System.out.println("2: liml,limu="+liml+" "+limu);        		 			
			liml=cbn1[je];
			limu=cbn2[je];
			//System.out.println("3: liml,limu="+liml+" "+limu);
			//System.out.println("ix1,ix2,ix3:"+ix1+" "+ix2+" "+ix3);
        		if((npr+limu-liml+1)>nprt){
       				nprt=nprt+5000;
				pls00=pls;
				plsnn=Array.getLength(pls00);
       				pls= new int[2*nprt];
				System.arraycopy(pls00,0,pls,0,plsnn);
				pls00=null;
        		}

			for( jj=liml;jj<=limu;jj++) {
				j=cbal[jj-1];
				radj=r0[j-1];
				if(radj>radprb&&j>i){
					ctf=rad+radj;
					ctf2=ctf*ctf;
					dctf=(float)Math.abs(rad-radj);
					dctf2=dctf*dctf;
					je=newIndexTwo(j,1,3);
					dx1=crd[je]-x1;
					dx2=crd[je+1]-x2;
					dx3=crd[je+2]-x3;
					d2=dx1*dx1+dx2*dx2+dx3*dx3;
 
	 				del=ctf2-d2;
	 				if(del>0.01&&d2>dctf2){
						npr=npr+1;
						je=newIndexTwo(npr,1,2);
						pls[je]=i;
						pls[je+1]=j;
						//System.out.println("pls "+npr+" "+pls[je]+" "+pls[je+1]);
					}
				}
			}
			if(npr==nprp){
				ast[i-1]=0;
			}
			nprp=npr;
		}
	 	System.out.println("npr = "+npr);

 
		cbln=rdmx+radprb;
		cubedata(2.0f,natm,cbln);
        	cbn1= new int[(lcb+1)*(mcb+1)*(ncb+1)];
        	cbn2= new int[(lcb+1)*(mcb+1)*(ncb+1)];
        	cbal= new int[27*natm];
		for(ip=1;ip<=1;ip++) {
			//je=newIndexTwo(ip,1,2);
			//System.out.println("pls2 "+ip+" "+npr+" "+pls[je]+" "+pls[je+1]);
		}
		cube(natm,crd,r0,radprb);

		for(int ix=0;ix<lcb+1;ix++) 
		for(int iy=0;iy<mcb+1;iy++) 
		for(int iz=0;iz<ncb+1;iz++) 
		{
			//je=ix*(mcb+1)*(ncb+1)+iy*(ncb+1)+iz;
			//System.out.println("cbn1,cbn2 -->"+ix+" "+iy+" "+iz+" "+(cbn1[je])+" "+cbn2[je]);
		}
		/*for(int ix=0;ix<lcb+1;ix++) 
		for(int iy=0;iy<mcb+1;iy++) 
		for(int iz=0;iz<ncb+1;iz++) 
		{
			je=ix*(mcb+1)*(ncb+1)+iy*(ncb+1)+ix;
			System.out.println("cbn1,cbn2 "+ix+" "+iy+" "+iz+" "+(cbn1[je])+" "+cbn2[je]);
		}*/	
		for(ip=1;ip<=1;ip++) {
			//je=newIndexTwo(ip,1,2);
			//System.out.println("pls1 "+ip+" "+npr+" "+pls[je]+" "+pls[je+1]);
		}
		nprx=0;
		float expos00[];
		int   exposnn=0;
		for(ip=1;ip<=npr;ip++) {
			je=newIndexTwo(ip,1,2);
			i=pls[je];
			j=pls[je+1];
			je=newIndexTwo(j,1,3);
			ke=newIndexTwo(i,1,3);
			//System.out.println("ip="+ip+" "+i+" "+j+" "+je+" "+ke);
			dx1=crd[je]-crd[ke];
			dx2=crd[je+1]-crd[ke+1];
			dx3=crd[je+2]-crd[ke+2];
			d2=dx1*dx1+dx2*dx2+dx3*dx3;
			dmg=(float)Math.sqrt(d2);
			pre=1.f+(r02[i-1]-r02[j-1])/d2;
			ke=newIndexTwo(i,1,3);
			tij1=crd[ke]+0.5f*pre*dx1;
			tij2=crd[ke+1]+0.5f*pre*dx2;
			tij3=crd[ke+2]+0.5f*pre*dx3;
			rij=0.5f*(float)Math.sqrt((r0[i-1]+r0[j-1])*(r0[i-1]+r0[j-1])-d2)*(float)
					Math.sqrt(d2-(r0[i-1]-r0[j-1])*(r0[i-1]-r0[j-1]))/dmg;
			//System.out.println("i,j,r0[i-1],r0[j-1],d2,dmg "+i+" "+j+" "+(r0[i-1])+" "+(r0[j-1])
			//	+" "+d2+" "+dmg+" "+rij);
 			//System.out.println("tij1,tij2,tij3,rij "+tij1+" "+tij2+" "+tij3+" "+rij);
			rvmg=(float)Math.sqrt(dx1*dx1+dx2*dx2);
			if(rvmg>1.0e-8){
				rv1=-dx2/rvmg;
				rv2=dx1/rvmg;
				cst=dx3/dmg;
				snt=(float) Math.sqrt(1.-cst*cst);
 
				csp=1.0f-cst;
				tm=csp*rv1;
				sm1=snt*rv1;
				sm2=snt*rv2;

				rm[1-1]=tm*rv1+cst;
				rm[4-1]=tm*rv2;
				rm[7-1]=sm2;
				rm[2-1]=tm*rv2;
				rm[5-1]=csp*rv2*rv2+cst;
				rm[8-1]=-sm1;
				rm[3-1]=-sm2;
				rm[6-1]=sm1;
				rm[9-1]=cst;
			}else {
        			rm[1-1]=1.0f;
        			rm[4-1]=0.0f;
        			rm[7-1]=0.0f;
       				rm[2-1]=0.0f;
        			rm[5-1]=1.0f;
        			rm[8-1]=0.0f;
        			rm[3-1]=0.0f;
        			rm[6-1]=0.0f;
        			rm[9-1]=1.0f;
			}
			nvo=0;
			nbv=0;
			// assign memory to expos if needed
        		if((nacc+nv)>nacct){
        			nacct=nacct+1000;
				expos00=expos;
        			expos=new float[3*nacct];
				if(expos00!=null) {
					exposnn=Array.getLength(expos00);
					System.arraycopy(expos00,0,expos,0,exposnn);
					expos00=null;	
				}			
        		}
			re10:
			for( iv=1;iv<=nvi;iv++) {
				je=newIndexTwo(iv,1,3);
				cf1=rm[1-1]*ver[je]+rm[4-1]*ver[je+1]+rm[7-1]*ver[je+2];
				cf2=rm[2-1]*ver[je]+rm[5-1]*ver[je+1]+rm[8-1]*ver[je+2];
				cf3=rm[3-1]*ver[je]+rm[6-1]*ver[je+1]+rm[9-1]*ver[je+2];
				cf1=tij1+rij*cf1;
				cf2=tij2+rij*cf2;
				cf3=tij3+rij*cf3;
 				
				ic1=(int)((cf1-xo)*cbai);
				ic2=(int)((cf2-yo)*cbai);
				ic3=(int)((cf3-zo)*cbai);
 				je=ic1*(mcb+1)*(ncb+1)+ic2*(ncb+1)+ic3;
				//System.out.println("cf1,cf2,cf3 "+iv+" "+cbai+" "+cf1+" "+cf2+" "+cf3+" "+xo+" "+yo+" "+zo);
        			//liml=cbn1[ic1+1+(lcb+1)*ic2+(lcb+1)*(mcb+1)*ic3-1];
        			//limu=cbn2[ic1+1+(lcb+1)*ic2+(lcb+1)*(mcb+1)*ic3-1];
				liml=cbn1[je];
				limu=cbn2[je];
				//System.out.println("liml,limu "+ic1+" "+ic2+" "+ic3+" "+liml+" "+limu);
				for(ii=liml;ii<=limu;ii++){
					k=cbal[ii-1];
					je=newIndexTwo(k,1,3);
					dy1=crd[je]-cf1;
					dy2=crd[je+1]-cf2;
					dy3=crd[je+2]-cf3;
					ds2=dy1*dy1+dy2*dy2+dy3*dy3;
					if(ds2<rs2[k-1]){
						//System.out.println("iv,k "+iv+" "+k);
						oti[iv-1]=k;
						continue re10;
					}
				}
				nvo=nvo+1;
				nacc=nacc+1;
 
				je=newIndexTwo(nacc,1,3);
				expos[je]=cf1;
				expos[je+1]=cf2;
				expos[je+2]=cf3;
				//System.out.println("nacc "+nacc+" "+expos[je]+" "+expos[je+1]+" "+expos[je+2]);
				oti[iv-1]=0;
			}

			nst=0;
			if(nlvl>0){
				for( ie=nvi;ie>=1;ie=ie-1) {
					je=newIndexTwo(ie,1,2);
					ia1=oti[edgv[je]-1];
					ia2=oti[edgv[je+1]-1];
					//int jef=edgv[je];
					//int jep=edgv[je+1];
					//System.out.println("ie,ia1,ia2:"+ie+" "+ia1+" "+ia2+" "+jef+" "+jep);
					if(ia1>0&&ia1==ia2)continue;
					nst=nst+1;
					st[nst-1]=ie;
					//System.out.println("nst, ie "+ie+" "+nst);
				}
			}
			
			if(nst>0){
				//**30	
				re30:
				while(true) {
					ie=st[nst-1];
					//System.out.println("nst, ie"+nst+" "+ie);
					nst=nst-1;
					je=newIndexTwo(ie,1,2);
					//System.out.println("nst, je"+nst+" "+je);
					ia1=oti[edgv[je]-1];
					ia2=oti[edgv[je+1]-1];
					iv=ie+nvi;
					je=newIndexTwo(iv,1,3);
					cf1=rm[1-1]*ver[je]+rm[4-1]*ver[je+1]+rm[7-1]*ver[je+2];
					cf2=rm[2-1]*ver[je]+rm[5-1]*ver[je+1]+rm[8-1]*ver[je+2];
					cf3=rm[3-1]*ver[je]+rm[6-1]*ver[je+1]+rm[9-1]*ver[je+2];
					cf1=tij1+rij*cf1;
					cf2=tij2+rij*cf2;
					cf3=tij3+rij*cf3;
 
					//if(ia1==0)goto 40;
					if(!(ia1==0)) {
						je=newIndexTwo(ia1,1,3);
						dy1=crd[je]-cf1;
						dy2=crd[je+1]-cf2;
						dy3=crd[je+2]-cf3;
						ds2=dy1*dy1+dy2*dy2+dy3*dy3;
 
						if(ds2<rs2[ia1-1]){
							oti[iv-1]=ia1;
							if(edg[ie-1]>0){
								nst=nst+1;
								st[nst-1]=edg[ie-1]+1;
							}
							//goto 60;
							if(nst>0) continue re30;
							else	  break re30;
						}
					}
					//40	continue
					//if(ia2==0)goto 50;
					if(!(ia2==0)) {
						je=newIndexTwo(ia2,1,3);
						dy1=crd[je]-cf1;
						dy2=crd[je+1]-cf2;
						dy3=crd[je+2]-cf3;
						ds2=dy1*dy1+dy2*dy2+dy3*dy3;
 
						if(ds2<rs2[ia2-1]){
							oti[iv-1]=ia2;
							if(edg[ie-1]>0){
								nst=nst+1;
								st[nst-1]=edg[ie-1];
							}
							//goto 60;
							if(nst>0) continue re30;
							else	  break re30;
						}
					}
					//50	continue
					ic1=(int)((cf1-xo)*cbai);
					ic2=(int)((cf2-yo)*cbai);
					ic3=(int)((cf3-zo)*cbai);
 
        				//liml=cbn1[ic1+1+(lcb+1)*ic2+(lcb+1)*(mcb+1)*ic3-1];
        				//limu=cbn2[ic1+1+(lcb+1)*ic2+(lcb+1)*(mcb+1)*ic3-1];
					je=ic1*(mcb+1)*(ncb+1)+ic2*(ncb+1)+ic3;
					liml=cbn1[je];
					limu=cbn2[je];
					for( ii=liml;ii<=limu;ii++) {
						k=cbal[ii-1];
						je=newIndexTwo(k,1,3);
						dy1=crd[je]-cf1;
						dy2=crd[je+1]-cf2;
						dy3=crd[je+2]-cf3;
						ds2=dy1*dy1+dy2*dy2+dy3*dy3;
						if(ds2<rs2[k-1]){
 
							oti[iv-1]=k;
							if(edg[ie-1]>0){
								nst=nst+1;
								st[nst-1]=edg[ie-1]+1;
								nst=nst+1;
								st[nst-1]=edg[ie-1];
							}
							if(nst>0) continue re30;
							else	  break re30;
							//goto 60;
						}
					}
					nvo=nvo+1;
					nacc=nacc+1;
 
					je=newIndexTwo(nacc,1,3);
					expos[je]=cf1;
					expos[je+1]=cf2;
					expos[je+2]=cf3;
					oti[iv-1]=0;
					if(edg[ie-1]>0){
						if(edg[edg[ie-1]+1-1]>0||ia2>0){
							nst=nst+1;
							st[nst-1]=edg[ie-1]+1;
						}
						if(edg[edg[ie-1]-1]>0||ia1>0){
							nst=nst+1;
							st[nst-1]=edg[ie-1];
						}
					}
				
					//60	if(nst>0)goto 30;
				
					if(nst>0) continue re30;
					else	  break re30;
				}
			}

			if(nvo>0){
				nprx=nprx+1;
				ast[i-1]=0;
				ast[j-1]=0;
			}

		}

		pls= null;
		cbn1= null;
		cbn2= null;
		cbal= null;

 

		nxa=0;
		for( i=1;i<=natm;i++) {
			if(ast[i-1]==0)nxa=nxa+1;
		}

        	System.out.println("no. pairs analyzed = "+npr);
        	System.out.println("no. exposed pairs  = "+nprx);
        	System.out.println("no. arc points  = "+nacc);
        	System.out.println("no. surface atoms  = "+nxa+" nbur = "+(natm-nxa));
		extot=nacc;
		for(i=1;i<=extot;i++) {
			je=newIndexTwo(i,1,3);
			//System.out.println(je);
 			//System.out.println("myexpose:"+i+" "+(expos[je])+" "+(expos[je+1])+" "+(expos[je+2]));
		}
	}





	public void qdiff () {
 
     
		String day;
	 	float  spec=0;	
 		int    noit;
		int    irm;
		long start;
		float  dbval[]=new float[2*7*2];
		float  sfd[]=new float[5*2];
		finish = cputime(0); 
		start=finish;
		wrt(1);
  
		
 
		// read in pdbfile and, if necesssary, assign charges/radii
		System.out.println("assigning atom information to delphi..."); 
		getatm();
	 	System.out.println("number of atoms for delphi.."+natom);
 		
		//finish=cputime(finish); 
		//System.out.println("time to read in and/or assign rad/chrg="+finish);
		 
		// find extrema
	 
		extrm();
 
		// atpos=atom positions, oldmid=current midpoints, rmaxdim=largest radius
		// rad3= radii in angstroms. calculate scale according to them and
		// to the percent box fill
 
		if(igrid==0){

			if(scale==10000.f)scale=1.2f;
			if(perfil==10000.)perfil=80.f;
			igrid=(int) (scale*100.f/perfil*rmaxdim);

		}else if(scale==10000.f){

			if(perfil==10000.f){
				scale=1.2f;
				perfil=100.f*rmaxdim*scale/(float)(igrid-1);
			}else {
				scale=(float)(igrid-1)*perfil/(100.f*rmaxdim);
			}

		}else {

			perfil=100.f*rmaxdim*scale/(float)(igrid-1);

		}
	 
		if(rionst>0.0f&&exrad<1.e-6f)exrad=2.0f;
	 
		irm=mod(igrid,2);
		if(irm==0)igrid=igrid+1;
		//
		if(igrid>ngrid){
			System.out.println("igrid = "+igrid+" exceeds maxgrid = "+ngrid) ;
     			System.out.println("so reset");
			igrid=ngrid;
			scale=(float)(igrid-1)*perfil/(100.f*rmaxdim);
		}
	 
		wrtprm();
	 
		ngp=igrid*igrid*igrid+1;
		nhgp=ngp/2;
		nbgp=igrid*igrid+1;
	
		iepsmp=new int[3*igrid*igrid*igrid];
		idebmap=new int[igrid*igrid*igrid];
 
		// calculate offsets
 
		off();
 
        	float xl1=oldmid[0]-(1.0f/scale)*(float)(igrid+1)*0.5f;
        	float xl2=oldmid[1]-(1.0f/scale)*(float)(igrid+1)*0.5f;
        	float xl3=oldmid[2]-(1.0f/scale)*(float)(igrid+1)*0.5f;
        	float xr1=oldmid[0]+(1.0f/scale)*(float)(igrid+1)*0.5f;
        	float xr2=oldmid[1]+(1.0f/scale)*(float)(igrid+1)*0.5f;
        	float xr3=oldmid[2]+(1.0f/scale)*(float)(igrid+1)*0.5f;
	 
		if(cmin[0]<xl1||cmin[1]<xl2||cmin[2]<xl3|| 
		   cmax[0]>xr1||cmax[1]>xr2||cmax[2]>xr3) {
    		 	System.out.println("!!! WARNING: part of molecule outside the box; hope that is OK");
		}

		// convert atom coordinates from angstroms to grid units
 	
		xn2=new float[3*natom];
		atmcrg=new float[4*natom];
		chgpos=new float[3*natom];
		crgatn=new int[natom]; 
		grdatm();

 
		// make some charge arrays for boundary conditions etc.
 
		if(isolv){
 
     			crgarr();
//
			if(logs){
				int ico=0;
				for( int ic=1;ic<=nqass;ic++) {
					int je=newIndexTwo(ic,1,3);
					float cx1=chgpos[je];
					float cx2=chgpos[je+1];
					float cx3=chgpos[je+2];
					if(cx1<xl1||cx1>xr1||cx2<xl2|| cx2>xr2||cx3<xl3||cx3>xr3){
						System.out.println("!!! FATAL ERROR: charge  outside the box");     
						ico=1;
					}
				}
				if(ico>0)return;
			}
		}
	 
		// write details
		wrtadt();
 
 
		// make the epsmap, and also a listing of boundary elements, and
		// the second epsmap used for the molecular surface scaling
 
		epsmak();
	
 
		if(isolv){
     			System.out.println("number of dielectric boundary points"+ibnum);
			if(iexun&&(ibnum==0)) {
				System.out.println("exiting as no boundary elements and");
				System.out.println("uniform dielectric exit flag has been set");
				return;
			}
		}
		//###################
		for(int i=0;i<=1;i++) 
		for(int j=0;j<=6;j++) 
		for(int k=0;k<=1;k++)
		{
			int je=i*7*2+j*2+k;
			//System.out.println("dbval i,j,k "+i+" "+j+" "+k+" "+dbval[je]);
			dbval[je]=0;
		}
		dbsfd(dbval,sfd);
		for(int i=1;i<=5;i++) 		
		for(int k=0;k<=1;k++)
		{
			int je=(i-1)*2+k;
			System.out.println("dbval i,j,k "+i+" "+k+" "+sfd[je]);
		}
	  
        	nsp=ibnum+1000;
        	db= new float[6*nsp];
       	 	idpos=new int[nsp];
        	sf1= new float[nhgp];
        	sf2= new float[nhgp];

		//mkdbsf(ibnum,nsp,dbval,icount2a,icount2b,sfd);
		mkdbsf(dbval,sfd);
	 	 
 
		// make qval and other linear charge arrays for the solver
 
		phimap=new float[igrid*igrid*igrid];
        	if(isph) {
			//setfcrg(nqgrd,icount1a,icount1b);
			//setfcrg();
		}else {
			//setcrg(nqgrd,icount1a,icount1b);
			
			setcrg();
			System.out.println("icount1a "+icount1a+" "+icount1b);
       		} 
 
		finish=cputime(finish);
		System.out.println("iepsmp to db, and charging done at "+ finish);
		System.out.println("number of grid points assigned charge "+ icount1b);
 
		//write dielectric map
 
		/*if(epswrt){
        		imaxwrd = igrid/16 + 1;
			neps=new int[3*imaxwrd*igrid*igrid];
			keps= new int[imaxwrd*igrid*igrid];
			wrteps(imaxwrd);
			neps= null;
			keps= null;
		}*/ 
	 
 
 
		// calculate boundary conditions
	 
		//setbc(qplus,qmin,cqplus,cqmin,nqass);
	 	setbc();
 
        	phimap1=new float[nhgp];
        	phimap2= new float[nhgp];
        	phimap3= new float[ngp];
 
        	bndx= new float[nbgp];
        	bndy= new float[nbgp];
        	bndz= new float[nbgp];

		if(iuspec) {
			spec=uspec;
			System.out.println("using entered value for relaxation of: "+spec);
		}else {
       	 		spec=relfac(spec);
		} 
 
		noit=(int)(7.8f/Math.log(1.0 + Math.sqrt(1-spec)));
		System.out.println("estimated iterations to convergence "+noit);
		if(iautocon) nlit = noit;
 
        	bndx1= new float[nbgp];
        	bndx2= new float[nbgp];
        	bndx3= new float[nbgp];
        	bndx4= new float[nbgp];
 
		finish=cputime(finish);
		System.out.println("  ");
		System.out.println("setup time was (sec) "+finish);
		System.out.println("  ");
 
		// iterate
 		 
		if(nnit==0||rionst<1.e-6) {

			itit(spec);

		}
 
        	bndx1= null;
        	bndx2= null;
        	bndx3= null;
        	bndx4= null;

        	bndx= null;
        	bndy= null;
        	bndz= null;
 
        	phimap1= null;
        	phimap2= null;
        	phimap3= null;
 
		db= null;
		idpos= null;
		sf1= null;
		sf2= null;
 
		encalc();
	
		//if(isite) wrtsit(nqass,icount2b,atpos,chrgv4,natom,ibnum);
 
		finish=cputime(start);
 
		System.out.println("  ");
		System.out.println("total cpu time was (sec) "+finish);
		System.out.println("  ");
		day=getTimeStr();
		System.out.println("DelPhi exited at "+day);
	 
	}
 
	 
 

	int  fxhl(int ivbeg,int ivend,int itbeg,int itend,int ntot,int vtpnt[], 
		   int vtlen[], int vindx[],int itot,int  vtlst[],float vert[],
		   int vtot,int tmlst[],int imxvtx,int  imxtri, float vnor[]) 
	{

		//**implicit none;

		int intsiz;
		int mxvtx=150000, mxtri=300000;		
		int pe,ke,je;  
		int is, ir, i, j1, j2, j3, i1, i2, n, it, id, nsel[];
		int tv[]=new int[20],vl[]=new int[100];
		float a1,b1, c1, a2, b2, c2, disnor, xnor, ynor, znor, p1, p2, p3;
		int j,k,n1;

		ntot=itot;
		nsel=new int[ivend-ivbeg+1000];
		for(i=ivbeg;i<=ivend;i++) {
	   		 nsel[i-1]=0;
		}

		is=0;
		ir=0;
		for( i=ivbeg;i<=ivend;i++) {
	    		i1=vtpnt[i-1];
	    		i2=i1+vtlen[i-1]-1;

	    		n=0;
	    		for( j=i1;j<=i2;j++) {
				k=3*vtlst[j-1];
				j1=vindx[k-3];
				j2=vindx[k-2];
				j3=vindx[k-1];
				if(j1>i){
		    			nsel[j1-1]=nsel[j1-1]+1;
		   			n=n+1;
		    			vl[n-1]=j1;
				}
				if(j2>i){
		    			nsel[j2-1]=nsel[j2-1]+1;
		    			n=n+1;
		    			vl[n-1]=j2;
				}
				if(j3>i){
		    			nsel[j3-1]=nsel[j3-1]+1;
		    			n=n+1;
		    			vl[n-1]=j3;
				}
				// k= one of the triangles
	    		}

	    		it=0;
	   		for( j=1;j<=n;j++) {
				k=vl[j-1];
				if((nsel[k-1]!=0)&&(nsel[k-1]!=2)) {
		    			it=it+1;
		    			tv[it-1]=k;
				}
				nsel[k-1]=0;
	   		}

	    		if(it>0) is=is+1;

				// mend
	    		if(it==2) {
				n1=tv[1-1];
				id=1;
				for( j=i1;j<=i2;j++) {
		    			k=3*vtlst[j-1];
		    			j1=vindx[k-2];
		    			j2=vindx[k-2];
		     			j3=vindx[k-1];
		     			if(j1==n1) id=1;
		     			if(j3==n1) id=2;
		     			if(j2==n1) {
			 			if(j1==i) id=2;
		    			}else {
			 			if(j3==i) id=1;
		     			}
		 		}
	
		 		j1=i;
		 		j2=tv[1-1];
		 		j3=tv[2-1];
		 		je=newIndexTwo(j2,1,3);
		 		ke=newIndexTwo(j1,1,3);
		 		a1=vert[je] - vert[ke];
		 		b1=vert[je+1] - vert[ke+1];
		 		c1=vert[je+2]-vert[ke+2];
		 		je=newIndexTwo(j3,1,3);
		 		ke=newIndexTwo(j1,1,3);
		 		a2=vert[je]-vert[ke];
		 		b2=vert[je+1]-vert[ke+1];
		 		c2=vert[je+2]-vert[ke+2];

		 		xnor=c2*b1-c1*b2;
		 		ynor=a2*c1-a1*c2;
		 		znor=a1*b2-a2*b1;

		 		disnor=(float)Math.sqrt(xnor*xnor+ynor*ynor+znor*znor);
		 		xnor=xnor/disnor;
		 		ynor=ynor/disnor;
		 		znor=znor/disnor;

		 		je=newIndexTwo(j1,1,3);
		 		ke=newIndexTwo(j2,1,3);
		 		pe=newIndexTwo(j3,1,3);
		 		p1=xnor*vnor[je]+ynor*vnor[je+1]+znor*vnor[je+2];
		 		p2=xnor*vnor[ke]+ynor*vnor[ke+1]+znor*vnor[ke+2];
		 		p3=xnor*vnor[pe]+ynor*vnor[pe+1]+znor*vnor[pe+2];

		 		if((p1<0)&&(p2<0)&&(p3<0)) {
		     			id=2;
		 		}else {
		     			id=1;
		 		}

		 		itot=itot+1;
		 		if(id==1) {
		     			vindx[3*itot-3]=i;
		     			vindx[3*itot-2]=tv[1-1];
		     			vindx[3*itot-1]=tv[2-1];
		 		}
		 		if(id==2) {
		     			vindx[3*itot-3]=tv[2-1];
		     			vindx[3*itot-2]=tv[1-1];
		     			vindx[3*itot-1]=i;
		 		}
		 		ir=ir+1;
	     		}

	 	}

		ntot=itot-ntot;

		// remake vertex to triangles listing..
		mkvtl(ivbeg,ivend,itbeg,itot, vtpnt, vtlen, vindx, vtlst, tmlst, imxtri, imxvtx);
		this.itot=itot;
		
		return ntot;
	}
	  
	void dbsfd(float dbval[],float sfd[]){
 
		//dimension dbval[0:1,0:6,0:1],sfd[5,0:1];
 		float denom,difeps,sixeps,debfct,sixth;
		int je,iz,ix,iy;
		debfct = epsout/((deblen*scale)*(deblen*scale));
		difeps=epsin-epsout;
		sixeps=epsout*6.0f;
		sixth=1.0f/6.0f;
 
		if(rionst>0.) {
			for( iz=0;iz<=1;iz++) {
	    		for( ix=1;ix<=3;ix++) {
	    			denom= sixeps + ix*difeps + iz*debfct;
	    			je=0*7*2+ix*2+iz;
	    			dbval[je]= 0.0f;
	   			je=1*7*2+ix*2+iz;	
	    			dbval[je]= difeps/denom;
				je=(ix-1)*2+iz;
	    			sfd[je]=epsout/denom;
			}
			}
			for(iz=0;iz<=1;iz++) {
	    		for(ix=4;ix<=5;ix++) {
	    			denom= sixeps + ix*difeps + iz*debfct;
		 		je=1*7*2+ix*2+iz;	
	    			dbval[je]= 0.0f;
		 		je=0*7*2+ix*2+iz;
	    			dbval[je]= -difeps/denom;
				je=(ix-1)*2+iz;
	    			sfd[je]=epsin/denom;
			}
			}
		}else {
			for( iz=0;iz<=1;iz++) {
	    		for( ix=1;ix<=5;ix++) {
	    			denom= sixeps + ix*difeps ;
		 		je=0*7*2+ix*2+iz;	
	    			dbval[je]= (epsout/denom) -sixth ;
		 		je=1*7*2+ix*2+iz;	
	    			dbval[je]= (epsin/denom) - sixth;
			}
			}
		}
	} 


	
	void mkdbsf(float dbval[],float sfd[]){
 
	
			 
		int deb,it[]=new int[6],itemp,icgrid;
 		int je,ke,ieps,ibnum3,iv,ibnum2,ibnum1,iw,isgrid,idbs,i,j,k,ix,iy,iz;
		float sixeps,debfct,sixth,sixsalt,dbs,difeps;

		iepsv= new int[nsp];
		phimap3= new float[ngp];
 
 		debfct = epsout/((deblen*scale)*(deblen*scale));
        	difeps=epsin-epsout;
        	sixeps=epsout*6.0f;
		sixth=1.0f/6.0f;
		sixsalt=sixth*((1/(1+debfct/sixeps))-1.0f);
		isgrid=igrid*igrid;
 	
		ibnum1=0;
		ibnum2=nsp/2;
		idbs=0;
		dbs=0.f;
 
		/*if(idbwrt) {
			open[13,file=dbnam[:dblen]);
			write[13,*) "DELPHI DB FILE";
			write[13,*) "FORMAT NUMBER=1";
			write[13,*) "NUMBER OF BOUNDARY POINTS= ",ibnum;
		}*/
 
		for( ix=1;ix<=ibnum;ix++) {
	  		je=newIndexTwo(ix,1,3);
	  		i=ibgrd[je];
	  		j=ibgrd[je+1];
	  		k=ibgrd[je+2];
 
	  		it[1-1]=0;
	  		it[2-1]=0;
	  		it[3-1]=0;
	  		it[4-1]=0;
	  		it[5-1]=0;
	  		it[6-1]=0;
	  		ieps=0;
	 		je=newIndexFour(i,j,k,1,igrid);
	  		if(iepsmp[je]!=0) it[0]=1;
	  		if(iepsmp[je+1]!=0) it[1]=1;
	  		if(iepsmp[je+2]!=0) it[2]=1;
			je=newIndexFour(i-1,j,k,1,igrid);
	  		if(iepsmp[je]!=0) it[3]=1;
		 	je=newIndexFour(i,j-1,k,2,igrid);
	  		if(iepsmp[je]!=0) it[4]=1;
			je=newIndexFour(i,j,k-1,3,igrid);
	  		if(iepsmp[je]!=0) it[5]=1;
	  		ieps=it[1-1]+it[2-1]+it[3-1]+it[4-1]+it[5-1]+it[6-1];
 
			je=newIndexThree(i,j,k,igrid);
	  		deb=idebmap[je];
	  		if(deb==1) {
	  			idbs=idbs+1;
	  		}
 
	  		iw=isgrid*(k-1) + igrid*(j-1) + i;
	  		iv=(iw+1)/2;
	  		if(iw!=(2*iv)) {
				ibnum1=ibnum1+1;
				ibnum3=ibnum1;
	  		}else {
				ibnum2=ibnum2+1;
				ibnum3=ibnum2;
	  		}
 
	  		idpos[ibnum3-1] = iv;
	  		iepsv[ibnum3-1] = ieps;
			//0:1,0:6,0:1
			je=it[3]*2*7+ieps*2+deb;
			ke=newIndexTwo(ibnum3,1,6);
	  		db[ke]=dbval[je];
			je=it[0]*2*7+ieps*2+deb;
	  		db[ke+1]=dbval[je];
			je=it[4]*2*7+ieps*2+deb;
	  		db[ke+2]=dbval[je];
			je=it[1]*2*7+ieps*2+deb;
	  		db[ke+3]=dbval[je];
			je=it[5]*2*7+ieps*2+deb;
	  		db[ke+4]=dbval[je];
			je=it[2]*2*7+ieps*2+deb;
	  		db[ke+5]=dbval[je];
	 
 
			ke=newIndexTwo(ibnum3,1,6);
	  		dbs=dbs+db[ke];
		}
 
	 
		System.out.println("no. dielectric boundary points in salt = "+idbs);
 
		// floatign idpos and db,compressing to contingous space
 
		icount2a=ibnum1;
		icount2b=icount2a+ibnum2-(nsp/2);
		itemp=(nsp/2);float ddd=0;
		for( ix=icount2a+1;ix<=icount2b;ix++) {
			itemp=itemp+1;
			idpos[ix-1]=idpos[itemp-1];
			iepsv[ix-1]=iepsv[itemp-1];
			ke=newIndexTwo(ix,1,6);
			je=newIndexTwo(itemp,1,6);
			db[ke]=db[je];
			db[ke+1]=db[je+1];
			db[ke+2]=db[je+2];
			db[ke+3]=db[je+3];
			db[ke+4]=db[je+4];
			db[ke+5]=db[je+5];
			ddd+=db[ke]+db[ke+1]*2+db[ke+2]*3+db[ke+3]*4+db[ke+4]*5+db[ke+5]*6;
			//System.out.println("db[ke+5]  "+ix+" "+db[ke]+" "+db[ke+1]+" "+db[ke+2]+" "+db[ke+3]);
		}
 		System.out.println("dbbb " +ddd);
		// set saltmaps 1 and 2
		// NB phimap3 used as a dummy flat 65**3 array
 
		if(rionst>0.) {
			iw=1;
			for( iz=1;iz<=igrid;iz++) {
	  		for( iy=1;iy<=igrid;iy++) {
	   		for( ix=1;ix<=igrid;ix++){
		  		je=newIndexThree(ix,iy,iz,igrid);
		  		phimap3[iw-1]=sixth + idebmap[je]*sixsalt;
	      			iw=iw+1;
			}
			}
			}
 
			iy=0;
			icgrid=igrid*igrid*igrid;
			sf1[(icgrid+1)/2-1]=phimap3[icgrid-1];
 
			for( ix=1;ix<=icgrid-2;ix+=2) {
				iy=iy+1;
				sf1[iy-1]=phimap3[ix-1];
				sf2[iy-1]=phimap3[ix+1-1];
			}
 			float ttt=0;
        		for( ix=1;ix<=icount2a;ix++) {
				i=1;
				if(sf1[idpos[ix-1]-1]==sixth) i=0;
				sf1[idpos[ix-1]-1]=sfd[(iepsv[ix-1]-1)*2+i];
				ttt+=sf1[idpos[ix-1]-1]*ix;
			}
        		for( ix=icount2a+1;ix<=icount2b;ix++) {
				i=1;
				if(sf2[idpos[ix-1]-1]==sixth) i=0;
				sf2[idpos[ix-1]-1]=sfd[2*(iepsv[ix-1]-1)+i];
				ttt+=sf2[idpos[ix-1]-1]*ix;
			}
 			System.out.println("ttt+ "+ttt);
		}
 		
		iepsv= null;
		phimap3= null;
	}
	


	void  setfcrg(int icount1a,int icount1b){
	//void  setfcrg() {
 
		//	Mike Gilson's fancy charge distribution embedded into qdiffx.
		//	A combination of Kim Sharp's chrgit1 and Anthony Nicholls' chrgup.
		//       Grafting by Richard Friedman.
		//	Latest version:6/2/89. First Version:4/25/89.
 
		boolean isum[]=new boolean[3*ngrid],flag;		 	 
		float ddist2[]=new float[500], rsave[]=new float[500];
		float clist[]=new float[4*50*500];
		int   nlist[]=new int[500];
       	 	int   kb[]=new int[3],pe,kz,ky,kx;
		float dg[]=new float[3];
		int   gchrg2[]=new int[ngcrg];
		float crit = .0001f;
		
		float gchrgd[]=new float[ngcrg];
		float radmin=1.2f;
		float radmax=2.0f;
		float radstep = .005f;
 
		float sum,dist,dist2,zdist,ydist,rmax,rmax2,sixeps,debfct,Rmax;
		float fpoh,difeps,dmin,ddz,ddx,ddy,epsins6,ysum,zsum,xsum,Rmax2,xdist;
		int il,je,ilist,irmax,k,i,j,ir,ig,n,ke,iy,ix,iz,itemp2,itemp1,ib,i2,i1,itemp;
		int ichoice,iv,isgrid,iw,itemp3;
		boolean goto1000,goto1500;

		Rmax2=0;Rmax=0;ichoice=0;
		// assign fractional charges to grid points, using phimap
		// as a dummy array (set to zero again later)
		difeps=epsin-epsout;
		sixeps=epsout*6.0f;
		fpoh=fpi*scale;
		debfct = epsout/((deblen*scale)*(deblen*scale));
		gchrgp=new int [3*ngcrg];
        	for( ig=1;ig<=nqgrd;ig++) {
        		for( i = 1;i<=3;i++) {
  
				// trucate to nearest grid point
			 	je=newIndexTwo(ig,i,4);
         		 	kb[i-1] = (int) (chrgv2[je]);
 
				// find position of charge in box in terms
				// of fractional distance along edges
 
          			dg[i-1] = chrgv2[je] - kb[i-1];
        	 
			}
 
			// loop over increasing radius to find suitable one
 
	  		irmax = 0;goto1000=false;
	  		for( rmax = radmin; rmax<=radmax; rmax+=radstep) {
	    			rmax2 = rmax*rmax;;
	    			irmax = irmax + 1;
	    			rsave [irmax-1] = rmax;
	    			ilist = 0;
	    			ir = (int)(rmax + dg[1-1] + 1);
 
				// loop over all grid points within range of rmax:
				// place base point of grid box containing the charge, at the origin, for now.
 
	    			for( i = -ir; i<=ir;i++) {
	      				xdist = i - dg[1-1];
	     				for( j = -ir; j<=ir;j++) {
	        				ydist = j - dg[2-1];
	        				for( k = -ir; k<=ir; k++) {
	         					zdist = k - dg[3-1];
 
							// calculate distance from current grid point to charge:
 
	          					dist2 =  xdist*xdist + ydist*ydist + zdist*zdist;
	          					dist = (float)Math.sqrt(dist2);
 
							// if grid point is closer than R, store the point and its distance:
 
		    					if (dist2<= Rmax2) {
	     							ilist = ilist + 1;
								je=(ilist-1)*500*4+(irmax-1)*4+1-1;
	     							clist[je]  = i;
	     							clist[je+1] = j;
	     							clist[je+2] = k;
	     							clist[je+3] = dist;
	          					}

						}
					}
				}
	   			nlist[irmax-1] = ilist;
 
				// generate weighting factors for this rmax:
				// sum normalizes the weighting to one:
 
	    			sum = 0;
	    			for(il = 1; il<=nlist[irmax-1];il++) {
					je=(il-1)*500*4+(irmax-1)*4+4-1;
	      				clist[je] = Rmax - clist[je];
	      				sum = sum + clist[je];
 
				}
 
				// normalize the weighting factors to sum to one:
 
	    			for( il = 1; il<=nlist[irmax-1];il++) {
					je=(il-1)*500*4+(irmax-1)*4+4-1;
	      				clist[je] = clist[je]/sum;
 
				}
 
				// calculate center of charge for this rmax:
 
	    			xsum = 0;
	    			ysum = 0;
	    			zsum = 0;
	    			for(il = 1; il<=nlist[irmax-1];il++) {
					je=(il-1)*500*4+(irmax-1)*4+1-1;
	      				xsum = xsum + clist[je+3] * clist[je];
	     		 		ysum = ysum + clist[je+3] * clist[je+1];
	      				zsum = zsum + clist[je+3] * clist[je+2];
 
				}
 
				// check whether criterion is satisfied, and if so, exit rmax loop.
 
	    			ddx = dg[1-1] - xsum;
	    			ddy = dg[2-1] - ysum;
	    			ddz = dg[3-1] - zsum;
	    			ddist2[irmax-1] = ddx * ddx + ddy*ddy + ddz * ddz;
	    			if (ddist2[irmax-1] <=crit) {
					
					goto1000=true;
					break;
 				}
				// otherwise, try another cutoff radius:		
 
 
			}
 
			// if loop gets finished without a radius yielding good enough
			// results, print warning,and use the best cutoff radius found:
			//	print *, 'Criterion not satisfied for charge #, ', iq
			//	print *, ' Using best cutoff radius found...'
 			goto1500=false;
 			if(!goto1000) {
 

				dmin = 100000;
				for(i = 1; i<=irmax;i++) {
	  				if (ddist2[i-1]<dmin) {
	    					dmin = ddist2[i-1];
	    					ichoice = i;
	  				}
				}

 
	  			dmin = (float) Math.sqrt (dmin);
	  			rmax = rsave[ichoice];
			
	  			goto1500=true;

			 
			}
			if(!goto1500) {
	  			ichoice = irmax;
	  			dmin = (float) Math.sqrt(ddist2[ichoice-1]);
			}
			//**1500	  continue

 
			// now we know what set of grids we're distributing the charge over
			// (clist(1-3, 1-ilist(ichoice), ichoice) for ichoice), and the weighting
			// (clist(4,1-ilist, ichoice)...
			// now, distribute the charge:
		 
	  		for( ilist = 1;  ilist<=nlist[ichoice-1]; ilist++) {
 
				// get grid point by adding offset to it
 		
				je=(ilist-1)*500*4+(ichoice-1)*4+1-1;
	    			kx =(int) ( kb[1-1] + clist[je]);
	    			ky =(int) ( kb[2-1] + clist[je+1]);
	    			kz =(int) (kb[3-1] + clist[je+2]);
 
				// make sure grid point is within the big box (should be a problem only
				// in cases where box edge cuts through
				// or very near the protein):
 
	    			flag = false;
   	    			if (kx<1 ) {
	      				kx = 1;
	      				flag = true;
	    			}
   	    			if (ky <1) {
	      				ky = 1;
	      				flag = true;
	    			}
   	    			if (kz <1) {
	      				kz = 1;
	      				flag = true;
	    			}
   	    			if (kx>igrid) {
	      				kx = igrid;
	      				flag = true;
	    			}
   	    			if (ky>igrid){
	      				ky = igrid;
	      				flag = true;
	    			}
   	    			if (kz>igrid) {
	      				kz = igrid;
	      				flag = true;
	    			}
		
	    			if(flag){
					
					for(j=1;j<=3;j++) {
						pe=newIndexTwo(ig,j,4);	
	      					System.out.println(" problem for charge at"+ chrgv2[pe]);
					}
	    			}
				je=newIndexThree(kx,ky,kz,igrid);
				ke=(ilist-1)*500*4+(ichoice-1)*4+4-1;pe=newIndexTwo(ig,4,4);
	   			 phimap[je]= phimap[je] + clist[ke]* chrgv2[pe];

 
			}

		}
		// set up odd/even boolean array
		for( i=1;i<=3*igrid;i+=2) {
			isum[i-i]=true;
		}
		for( i=2;i<=(3*igrid-1);i+=2) {
 			isum[i-i]=false;
		}
 
		// find which grid points have charge assigned to them
		// (will use this array later to calculate grid energy)
 
		n=0;
        	for( k=2;k<=igrid-1;k++) {
          	for( j=2;j<=igrid-1;j++){
            	for(i=2;i<=igrid-1;i++) {
			je=newIndexThree(i,j,k,igrid);
                	if(phimap[je]!=0)n=n+1;
	    	}
	  	}
		}
        	gchrgp=new int[3*n];
        	gchrg=new float[n];
        	gchrgd=new float[n];
        	gchrg2=new int[n];
		qval=new float[n];
        	iqpos=new int[n];

		for( k=2;k<=igrid-1;k++) {
	  	for( j=2;j<=igrid-1;j++) {
	    	for( i=2;i<=igrid-1;i++) {
			ke=newIndexThree(i,j,k,igrid);
			if(phimap[ke]!=0) {
				n=n+1;
				je=newIndexTwo(n,1,3);
				gchrgp[je]=i;
				gchrgp[je+1]=j;
				gchrgp[je+2]=k;
				gchrg[n-1]=phimap[ke];
				phimap[ke]=0.0f;
			}
		}
		}
		}

		icount1b=n;
 
		// determine how many charged grid points are odd
 
		icount1a=0;
		for( i=1;i<=n;i++) {
			je=newIndexTwo(i,1,3);
			itemp=gchrgp[je]+gchrgp[je+1]+gchrgp[je+2];
			if(isum[itemp-1]) icount1a=icount1a+1;
		}
 
		// set up odd/even pointer array, to be used in making qval
		// and iqpos
 
		i1=0;
		i2=icount1a;
		for( i=1;i<=n;i++) {
			je=newIndexTwo(i,1,3);
			itemp=gchrgp[je]+gchrgp[je+1]+gchrgp[je+2];
			if(isum[itemp-1]) {
				i1=i1+1;
				gchrg2[i-1]=i1;
			}else {
				i2=i2+1;
				gchrg2[i-1]=i2;
			}
		}
 
		// determine denominator at all charged grid points
 
		ib=0;
		epsins6=sixeps+ 6.0f*difeps;
		for( i=1;i<=n;i++){
			je=newIndexTwo(i,1,3);
			iz=gchrgp[je];
			iy=gchrgp[je+1];
			ix=gchrgp[je+2];
			je=newIndexFour(ix,iy,iz,1,igrid);
			ke=newIndexFour(ix-1,iy,iz,1,igrid);
			itemp1=iepsmp[je]+iepsmp[ke];
			ke=newIndexFour(ix,iy-1,iz,2,igrid);
			itemp2=iepsmp[je+1]+iepsmp[ke];
			ke=newIndexFour(ix,iy,iz-1,3,igrid);
			itemp3=iepsmp[je+2]+iepsmp[ke];
			itemp=itemp1+itemp2+itemp3;
			je=newIndexThree(ix,iy,iz,igrid);
			gchrgd[i-1]=itemp*difeps + debfct*idebmap[je] + sixeps;
			if(itemp!=6) {
				ib=ib+1;
				je=newIndexTwo(ib,1,5);
				ke=newIndexThree(ix,iy,iz,igrid);
				cgbp[je]=ix;
				cgbp[je+1]=iy;
				cgbp[je+2]=iz;
				cgbp[je+3]=gchrg[i-1]*fpoh/(epsins6+debfct*idebmap[ke]);
				cgbp[je+4]=gchrg2[i-1];
			}
		}
		ibc=ib;
		System.out.println("# grid points charged and at boundary="+ib);
 
		// make qval, fpoh term so potentials will be in kt/e
 
		for( i=1;i<=n;i++){
			j=gchrg2[i-1];
			qval[j-1]=gchrg[i-1]*fpoh/gchrgd[i-1];
		}
 
		// make iqpos
 
		isgrid=igrid*igrid;
		for( i=1;i<=n;i++) {
			j=gchrg2[i-1];
			je=newIndexTwo(i,1,3);
			ix=gchrgp[je];
			iy=gchrgp[je+1];
			iz=gchrgp[je+2];
			iw=1+ix+igrid*(iy-1)+isgrid*(iz-1);
			iv=iw/2;
			iqpos[j-1]=iv;
		}
 
		// end of chrgup, return with qval,iqpos and gchrgp and gchrg
		// also icount1a, icount1b
	}


	void  setbc(){
 
		// assigns either zero boundary conditions (ibctyp=1)
		// quasi-coulombic based on debye dipole qplus at cqplus
		// and qmin at cqmin (ibctyp=2)
		// focussing (ibctyp=3)(continuance if scales are the same)
		// full quasi-coulombic (ibctyp=4)
		// or constant external filed (ibctyp=5)
		// option 2 will be appropriately modified by periodic
		// boundary condition flags (iper)
 
	 	int je,ix,iy,iz,midg,ic;
		float g[]=new float[3],go[]=new float[3],c[]=new float[3];
		String toblbl,label,title,filnam,botlbl;
 		float dist,tempd,tempp,tempn,dist2,dist1,dist4,dist3;
		// zero option, clear boundary values
 
		System.out.println(" ");
		System.out.println("  setting boundary conditions");
		System.out.println(" ");
 
          	for( iz=1;iz<=igrid;iz++) {
               	for( ix=1;ix<=igrid;ix+=igrid-1){
             	for( iy=1;iy<=igrid;iy++) {
		   	je=newIndexThree(ix,iy,iz,igrid);
                   	phimap[je] = 0.0f;
		}
		}
		}
          	for(iz=1;iz<=igrid;iz++) {
              	for(iy=1;iy<=igrid;iy+=igrid-1){
               	for(ix=1;ix<=igrid;ix++) {
		 	je=newIndexThree(ix,iy,iz,igrid);
                   	phimap[je] = 0.0f;
                   	 
		}
		}
		}
          	for( iz=1;iz<=igrid;iz+=igrid-1) {
                for( iy=1;iy<=igrid;iy++) {
                for( ix=1;ix<=igrid;ix++) {
                   	je=newIndexThree(ix,iy,iz,igrid);
                   	phimap[je] = 0.0f;
                   	
		}
		}
		}
 
		// end of zero option
 
      		if(ibctyp==5) {
 
			// constant field, 1 kT/e/grid unit, in the x direction
 
         		for( iz=1;iz<=igrid;iz++) {
           		for( iy=1;iy<=igrid;iy++) {
               		for( ix=1;ix<=igrid;ix++) {
		 		je=newIndexThree(ix,iy,iz,igrid);
                   		phimap[je] =ix;
                  	}}}
		}
 
		// end external field option
 
      		if(ibctyp==2) {
 
			// quasi coulombic dipole option
 
			for( iz=1;iz<=igrid;iz++){
	  		for( iy=1;iy<=igrid;iy++) {
            			dist1 = (cqplus[1-1]-1)*(cqplus[1-1]-1) 
				      + (cqplus[2-1]-iy)*(cqplus[2-1]-iy)
				      + (cqplus[3-1]-iz)*(cqplus[3-1]-iz);
 	      			dist1 = (float)Math.sqrt(dist1)/scale;
 	      			tempp = qplus*(float)Math.exp(-dist1/deblen )/(dist1*epsout) ;
            			dist2 = (cqmin[1-1]-1)*(cqmin[1-1]-1)
					 + (cqmin[2-1]-iy)*(cqmin[2-1]-iy)
					 + (cqmin[3-1]-iz)*(cqmin[3-1]-iz);
 	     		 	dist2 = (float)Math.sqrt(dist2)/scale;
 	      			tempn = qmin*(float)Math.exp(-dist2/deblen )/(dist2*epsout) ;
		 		je=newIndexThree(1,iy,iz,igrid);
            			phimap[je] =  tempp + tempn;
            			dist3 = (cqplus[1-1]-igrid)*(cqplus[1-1]-igrid) 
				      + (cqplus[2-1]-iy)*(cqplus[2-1]-iy)	
				      + (cqplus[3-1]-iz)*(cqplus[3-1]-iz);
 	      			dist3 = (float)Math.sqrt(dist3)/scale;
 	      			tempp = qplus*(float)Math.exp(-dist3/deblen )/(dist3*epsout);
            			dist4 = (cqmin[1-1]-igrid)*(cqmin[1-1]-igrid)
					 + (cqmin[2-1]-iy)* (cqmin[2-1]-iy)
						+ (cqmin[3-1]-iz)*(cqmin[3-1]-iz);
 	      			dist4 = (float)Math.sqrt(dist4)/scale;
 	      			tempn = qmin*(float)Math.exp(-dist4/deblen )/(dist4*epsout) ;
		 		je=newIndexThree(igrid,iy,iz,igrid);
            			phimap[je] =  tempp + tempn;
				//System.out.println("dist1 "+iz+" "+iy+" "+tempp+" "+dist2+" "+dist3+" "+dist4+" "+deblen+" "+epsout);
			}
			}
			for( iz=1;iz<=igrid;iz++) {
	  		for( iy=1;iy<=igrid;iy+=igrid-1) {
	   		for( ix=1;ix<=igrid;ix++) {
            			dist = (cqplus[1-1]-ix)*(cqplus[1-1]-ix) 
					+ (cqplus[2-1]-iy)*(cqplus[2-1]-iy)
					 + (cqplus[3-1]-iz)*(cqplus[3-1]-iz);
 	      			dist = (float)Math.sqrt(dist)/scale;
 	      			tempp = qplus*(float)Math.exp(-dist/deblen )/(dist*epsout) ;
            			dist = (cqmin[1-1]-ix)*(cqmin[1-1]-ix) + 
					(cqmin[2-1]-iy)*(cqmin[2-1]-iy) + (cqmin[3-1]-iz)*(cqmin[3-1]-iz);
 	      			dist = (float)Math.sqrt(dist)/scale;
 	      			tempn = qmin*(float)Math.exp(-dist/deblen )/(dist*epsout) ;
		 		je=newIndexThree(ix,iy,iz,igrid);
            			phimap[je] =  tempp + tempn;
			}
			}
			}
			for( iz=1;iz<=igrid;iz+=igrid-1) {
	  		for( iy=1;iy<=igrid;iy++) {
	    		for( ix=1;ix<=igrid;ix++) {
            			dist = (cqplus[1-1]-ix)*(cqplus[1-1]-ix)
					 + (cqplus[2-1]-iy)*(cqplus[2-1]-iy) + 
					  (cqplus[3-1]-iz)*(cqplus[3-1]-iz);
 	      			dist = (float)Math.sqrt(dist)/scale;
 	      			tempp = qplus*(float)Math.exp(-dist/deblen )/(dist*epsout) ;
            			dist = (cqmin[1-1]-ix)*(cqmin[1-1]-ix) 
					+ (cqmin[2-1]-iy)*(cqmin[2-1]-iy)
					 + (cqmin[3-1]-iz)*(cqmin[3-1]-iz);
 	      			dist = (float)Math.sqrt(dist)/scale;
 	      			tempn = qmin*(float)Math.exp(-dist/deblen )/(dist*epsout) ;
		 		je=newIndexThree(ix,iy,iz,igrid);
            			phimap[je] = tempp + tempn;
			}
			}
			}
 
      		} 
 
		// end of quasi coulombic dipole option
 
      		if(ibctyp==4) {
 
			// a summation of the potential resulted from each point of charge 
 
			for(iz=1;iz<=igrid;iz++) {
	  		for( iy=1;iy<=igrid;iy++) {
	    		for( ix=1;ix<=igrid;ix+=igrid-1){
			for( ic = 1;ic<=nqass;ic++) {
				je=newIndexTwo(ic,1,4);
              			dist = (atmcrg[je]-ix)*(atmcrg[je]-ix)
					 + (atmcrg[je+1]-iy)*(atmcrg[je+1]-iy)
					 + (atmcrg[je+2]-iz)*(atmcrg[je+2]-iz);
 	        		dist = (float)Math.sqrt(dist)/scale;
	        		tempd = atmcrg[je+3]*(float)Math.exp(-dist/deblen)/(dist*epsout);
				je=newIndexThree(ix,iy,iz,igrid);
              			phimap[je] = phimap[je] + tempd ;
			}}}}

			for(iz=1;iz<=igrid;iz++) {
	  		for( iy=1;iy<=igrid;iy+=igrid-1) {
	   		for(ix=1;ix<=igrid;ix++) {
			for(ic = 1;ic<=nqass;ic++) {
				je=newIndexTwo(ic,1,4);
              			dist = (atmcrg[je]-ix)*(atmcrg[je]-ix) + (atmcrg[je+1]-iy)*(atmcrg[je+1]-iy)
					 + (atmcrg[je+2]-iz)*(atmcrg[je+2]-iz);
 	        		dist = (float)Math.sqrt(dist)/scale;
	        		tempd = atmcrg[je+3]*(float)Math.exp(-dist/deblen)/(dist*epsout);
				je=newIndexThree(ix,iy,iz,igrid);
              			phimap[je] = phimap[je] + tempd ;
			}}}}
			for(iz=1;iz<=igrid;iz+=igrid-1) {
	  		for( iy=1;iy<=igrid;iy++) {
	      		for(ix=1;ix<=igrid;ix++) {
			for(ic = 1;ic<=nqass;ic++) {
				je=newIndexTwo(ic,1,4);
              			dist = (atmcrg[je]-ix)*(atmcrg[je]-ix) + (atmcrg[je+1]-iy)*(atmcrg[je+1]-iy)
					 + (atmcrg[je+2]-iz)*(atmcrg[je+2]-iz);
             	  		dist = (float)Math.sqrt(dist)/scale;
 	      	  		tempd = atmcrg[je+3]*(float)Math.exp(-dist/deblen)/(dist*epsout);
	       			je=newIndexThree(ix,iy,iz,igrid);
              			phimap[je] = phimap[je] + tempd ;
          	
			}}}}
      		} 
 
		// end of the option for the complete summation of potential
 
      		if(ibctyp==3) {
  
		} 
 
		// end of focussing option
		//
 
 
		midg = (igrid+1)/2;
		System.out.println(" some initial phi values: ");
		System.out.println(" midg,midg,1; midg,midg,igrid ");
		System.out.print((phimap[newIndexThree(midg,midg,1,igrid)])+" "+ (phimap[newIndexThree(midg,midg,igrid,igrid)]));
		System.out.print(" midg,1,midg; midg,igrid,midg ");
		System.out.println((phimap[newIndexThree(midg,1,midg,igrid)])+" "+(phimap[newIndexThree(midg,igrid,midg,igrid)]));
		System.out.print(" 1,midg,midg; igrid midg midg: ");
		je=newIndexThree(1,midg,midg,igrid);
		System.out.println((phimap[je])+"  "+(phimap[newIndexThree(igrid,midg,midg,igrid)]));
  		System.out.println("ibctyp "+ibctyp);
	}
 


	


	float relfac(float spec) {
 
	 
		float sn1[]=new float[ngrid],sn2[]=new float[ngrid],sn3[]=new float[ngrid];
		 
	 
		int star,fin,sta1[]=new int[ngrid],sta2[]=new int[ngrid];
		int fi1[]=new int[ngrid],fi2[]=new int[ngrid];
       		int ix,iy,iz,iw,itemp2,i,j,k,inc1zb,inc2zb,inc1za,idif1z,n,inc1yb,inc2ya,iadd2;
		int inc2xb,idif2y,itemp1,itemp3,itemp4,idif2z,ihgd2,je,long2,long1,lat2,lat1,iadd1; 
		int icgrid,inc2xa,inc2za,ihgd,inc1xb,idif1x,inc1ya,inc2yb,inc1xa,idif1y,idif2x;
	
		float temp2=0,temp3=0,temp1=0,temp4=0,temp,sixth;
 		inc2xb=0;idif2x=0;inc2yb=0;idif1x=0;idif2z=0;inc2xa=0;inc2ya=0;idif2y=0;inc2zb=0;
		inc2za=0;inc1xb=0;inc1xa=0;inc1yb=0;inc1ya=0;idif1y=0;inc1zb=0;inc1za=0;idif1z=0;
		sixth = 1.f/6.f;
		icgrid=igrid*igrid*igrid;
		ihgd=(igrid+1)/2;
 
		if(iper[0]) {
			n=0;
			for( iz=2;iz<=igrid-1;iz++) {
	  			iadd1=(iz-1)*igrid*igrid ;
	  			for( iy=2;iy<=igrid-1;iy++) {
	    				iadd2=(iadd1+(iy-1)*igrid +2)/2;
	    				n=n+1;
	   				bndx[n-1]=iadd2;
				}
			}
			idif1x=(igrid-2)/2;
			idif2x=idif1x+1;
			inc1xa=1;
			inc1xb=0;
			inc2xa=0;
			inc2xb=1;
		} 
		if(iper[1]) {
			n=0;
			for( iz=2;iz<=igrid-1;iz++) {
	  			iadd1=(iz-1)*igrid*igrid ;
	  			for(ix=2;ix<=igrid-1;ix++) {
	    				iadd2=(iadd1+ix+1)/2;
	    				n=n+1;
    	    				bndy[n-1]=iadd2;
				}
			}
			idif1y=igrid*(igrid-2)/2;
			idif2y=idif1y+1;
			inc1ya=(igrid/2)+1;
			inc1yb=inc1ya-1;
			inc2ya=inc1yb;
			inc2yb=inc1ya;
		}
		if(iper[2]) {
			n=0;
			for( ix=2;ix<=igrid-1;ix++) {
	  			iadd1=ix+1;
	  			for( iy=2;iy<=igrid-1;iy++) {
	    				iadd2=(iadd1+(iy-1)*igrid)/2;
	    				n=n+1;
	    				bndz[n-1]=iadd2;
				}
			}
			idif1z=igrid*igrid*(igrid-2)/2;
			idif2z=idif1z+1;
			inc1za=((igrid*igrid)/2)+1;
			inc1zb=inc1za;
			inc2za=inc1zb;
			inc2zb=inc1za;
		}
 
		// set up start and stop vectors
		sta1[2-1]=(igrid*igrid + igrid +4)/2;
		sta2[2-1]=sta1[2-1]-1;
		fi1[2-1]=igrid*igrid - (igrid+1)/2;
		fi2[2-1]=fi1[2-1];
		itemp1=igrid + 2;
		itemp2=igrid*igrid -igrid -2;
		for(i=3;i<=igrid-1;i++) {
			sta1[i-1]=fi1[i-1-1] + itemp1;
			sta2[i-1]=fi2[i-1-1] + itemp1;
			fi1[i-1]=sta1[i-1-1] + itemp2;
			fi2[i-1]=sta2[i-1-1] + itemp2;
		}
 
		// also
 
		lat1= (igrid-1)/2;
		lat2= (igrid+1)/2;
		long1= (igrid*igrid - 1)/2;
		long2= (igrid*igrid + 1)/2;
 
		// set up sn array for lowest eigenstate
 
      		i=0;
 
		sn1[1-1]=0.0f;
		sn1[igrid-1]=0.0f;
		sn2[1-1]=0.0f;
		sn2[igrid-1]=0.0f;
		sn3[1-1]=0.0f;
		sn3[igrid-1]=0.0f;
		for( ix=2;ix<=igrid-1;ix++) {
        		temp=3.14159265f*(float)(ix-1)/(float)(igrid-1);
        		sn1[ix-1]=(float)Math.sqrt(2.0)*(float)Math.sin(temp)/(float)Math.sqrt((float)(igrid-1));
			sn2[ix-1]=sn1[ix-1];
			sn3[ix-1]=sn1[ix-1];
		}
		if(iper[0]) {
			for( ix=1;ix<=igrid;ix++) {
				sn1[ix-1]=1.0f/(float)Math.sqrt((float)(igrid));
			}
		}
		if(iper[1]) {
			for( iy=1;iy<=igrid;iy++) {
 				sn2[ix-1]=1.0f/(float)Math.sqrt((float)(igrid));
			}
		}
		if(iper[2]) {
			for( iz=1;iz<=igrid;iz++) {
 				sn3[iz-1]=1.0f/(float)Math.sqrt((float)(igrid));
			}
		}
 
		iw=1;
		for( iz=1;iz<=igrid;iz++) {
	  		temp3=sn3[iz-1];
	  		for( iy=1;iy<=igrid;iy++){
	    			temp2=temp3*sn2[iy-1];
	   			for(ix=1;ix<=igrid;ix++) {
	       				phimap3[iw-1]=temp2*sn1[ix-1];
	    				iw=iw+1;
				}
			}
		}
		temp=0.0f;
		for( ix=2;ix<=icgrid-1;ix+=2) {
			iy=ix/2;
			phimap2[iy-1]=phimap3[ix-1];
			temp=temp + phimap3[ix-1]*phimap3[ix-1];
		}
 
		if(rionst>0.0) {
        		for( n = 2; n<=igrid-1;n++) {
				star=sta1[n-1];
				fin=fi1[n-1];
            			for( ix = star;ix<=fin;ix++) {
            				temp1 = phimap2[ix-1] +  phimap2[ix-1-1];
 
            				temp2 = phimap2[ix+lat1-1] +   phimap2[ix-lat2-1];
 
            				temp3 = phimap2[ix+long1-1] +    phimap2[ix-long2-1];
 
       					phimap1[ix-1] = (temp1+temp2+temp3)*sf1[ix-1];
				}
			}
 
			// otherwise the main loop is as below:
 
        	}else {
 
        		for( n = 2; n<=igrid-1;n++) {
	    			star=sta1[n-1];
	    			fin=fi1[n-1];
            			for( ix = star;ix<=fin;ix++) {
            				temp1 = phimap2[ix-1] +  phimap2[ix-1-1];
 
            				temp2 = phimap2[ix+lat1-1] +   phimap2[ix-lat2-1];
 
            				temp3 = phimap2[ix+long1-1] + phimap2[ix-long2-1];
 
       					phimap1[ix-1] =  (temp1+temp2+temp3)*sixth;
				}
			}
       	 	}
 
		//$DIR NO_RECURRENCE 
		for( k=1;k<=icount2a;k++) {
			ix=idpos[k-1];
			je=newIndexTwo(k,1,6);
			temp1=phimap2[ix-1-1]*db[je]+phimap2[ix-1]*db[je+1];
			temp2=phimap2[ix-lat2-1]*db[je+2]+phimap2[ix+lat1-1]*db[je+3];
			temp3=phimap2[ix-long2-1]*db[je+4]+phimap2[ix+long1-1]*db[je+5];
			phimap1[ix-1]= phimap1[ix-1] + temp1+temp2+temp3;
		}
 
		// Now reset boundary values altered in above loops.
 
		star=igrid*(igrid+1)/2;
		fin=igrid*(igrid*(igrid-1)-2)/2;
		//$DIR NO_RECURRENCE
		for( ix=star;ix<=fin;ix+=igrid) {
	  		phimap1[ix+1-1]=0.0f;
	  		phimap1[ix+ihgd-1]=0.0f;
		}
 
		temp=0.0f;
		for( ix=1;ix<=(icgrid-1)/2;ix++) {
			temp=temp + phimap1[ix-1]*phimap3[2*ix-1-1];
		}
		// if periodic boundary condition option
		// force periodicity using wrap around update of boundary values:
		// 2nd slice-->last
		// last-1 slice-->first
		//
		// z periodicity
		//
        	if(iper[2]) {
         		for( iz = 1;iz<=(igrid-1)*(igrid-1);iz+=2) { 
	  			temp1=bndz[iz-1];
	  			temp2=temp1+idif1z;
	  			temp3=temp2+inc1za;
	  			temp4=temp1+inc1zb;
        			itemp1=(int)temp1;
        			itemp2=(int)temp2;
        			itemp3=(int)temp3;
        			itemp4=(int)temp4;
            			phimap1[itemp1-1]=phimap2[itemp2-1];
            			phimap1[itemp3-1]=phimap2[itemp4-1];
			}
        	}
 
		// y periodicity
 
        	if(iper[1]) {
         		for( iy = 1;iy<=(igrid-1)*(igrid-1);iy+=2) {
	  			temp1=bndy[iy-1];
	  			temp2=temp1+idif1y;
	  			temp3=temp2+inc1ya;
	  			temp4=temp1+inc1yb;
        			itemp1=(int)temp1;
        			itemp2=(int)temp2;
        			itemp3=(int)temp3;
        			itemp4=(int)temp4;
           	 		phimap1[itemp1-1]=phimap2[itemp2-1];
            			phimap1[itemp3-1]=phimap2[itemp4-1];
			}
        	}
 
		// x periodicity
 
        	if(iper[0]) {
          		for( ix = 1;ix<=(igrid-1)*(igrid-1);ix+=2) {
	  			temp1=bndx[ix-1];
	  			temp2=temp1+idif1x;
	  			temp3=temp2+inc1xa;
	  			temp4=temp1+inc1xb;
        			itemp1=(int)temp1;
        			itemp2=(int)temp2;
        			itemp3=(int)temp3;
        			itemp4=(int)temp4;
            			phimap1[itemp1-1]=phimap2[itemp2-1];
            			phimap1[itemp3-1]=phimap2[itemp4-1];
			}
         	}
 
		// Next update phimap3 using the new phimap1
 
        	if(rionst>0.0) {	
        		for( n = 2; n<=igrid-1;n++) {
	    			star=sta2[n-1];
	    			fin=fi2[n-1];
           			for( ix = star;ix<=fin;ix++) {
            				temp1 = phimap1[ix-1] + phimap1[ix+1-1];
 
            				temp2 = phimap1[ix+lat2-1] +  phimap1[ix-lat1-1];
 
            				temp3 = phimap1[ix+long2-1] + phimap1[ix-long1-1];
 
       					phimap3[ix-1] =(temp1+temp2+temp3)*sf2[ix-1];
				}
			}
 
		}else {
 
       			for(n = 2; n<=igrid-1;n++){
	    			star=sta2[n-1];
	    			fin=fi2[n-1];
            			for( ix = star;ix<=fin;ix++) {
            				temp1 = phimap1[ix-1] + phimap1[ix+1-1];
 
            				temp2 = phimap1[ix+lat2-1] + phimap1[ix-lat1-1];
 
            				temp3 = phimap1[ix+long2-1] + phimap1[ix-long1-1];
 
       					phimap3[ix-1] = (temp1+temp2+temp3)*sixth;
				}
			}
		}
 
		//$DIR NO_RECURRENCE 
		for( k=icount2a+1;k<=icount2b;k++) {
			ix=idpos[k-1];
			je=newIndexTwo(k,1,6);
			temp1=phimap1[ix-1]*db[je]+phimap1[ix+1-1]*db[je+1];
			temp2=phimap1[ix-lat1-1]*db[je+2]+phimap1[ix+lat2-1]*db[je+3];
			temp3=phimap1[ix-long1-1]*db[je+4]+phimap1[ix+long2-1]*db[je+5];
			phimap3[ix-1]=phimap3[ix-1] + temp1+temp2+temp3;
		}
		// reset boundary condition
 
		star=(igrid+2)/2;
		iy=(igrid*(igrid+2)/2) - igrid +1;
		fin=(igrid*(igrid-1)-1)/2;
		ihgd2=ihgd-1;
		//$DIR NO_RECURRENCE
		for( ix=star;ix<=fin;ix++){
			iy=iy+igrid;
			phimap3[iy-1]=0.0f;
			phimap3[iy+ihgd2-1]=0.0f;
		}
 
		// z periodicity
 
        	if(iper[2]) {
          		for( iz = 2;iz<=(igrid-1)*(igrid-1);iz+=2) {
	  			temp1=bndz[iz-1];
	  			temp2=temp1+idif2z;
	  			temp3=temp2+inc2za;
	  			temp4=temp1+inc2zb;
				itemp1=(int)temp1;
				itemp2=(int)temp2;
				itemp3=(int)temp3;
				itemp4=(int)temp4;
	    			phimap3[itemp1-1]=phimap1[itemp2-1];
	    			phimap3[itemp3-1]=phimap1[itemp4-1];
			}
        	}
 
		// y periodicity
 
        	if(iper[1]) {
          		for( iy = 2;iy<=(igrid-1)*(igrid-1);iy+=2) {
	  			temp1=bndy[iy-1];
	  			temp2=temp1+idif2y;
	  			temp3=temp2+inc2ya;
	  			temp4=temp1+inc2yb;
        			itemp1=(int)temp1;
        			itemp2=(int)temp2;
        			itemp3=(int)temp3;
        			itemp4=(int)temp4;
            			phimap3[itemp1-1]=phimap1[itemp2-1];
            			phimap3[itemp3-1]=phimap1[itemp4-1];
			}
        	}
 
		// x periodicity
 
        	if(iper[0]) {
          		for( ix = 2;ix<=(igrid-1)*(igrid-1);ix+=2) {
	  			temp1=bndx[ix-1];
	  			temp2=temp1+idif2x;
	  			temp3=temp2+inc2xa;
	  			temp4=temp1+inc2xb;
        			itemp1=(int)temp1;
        			itemp2=(int)temp2;
        			itemp3=(int)temp3;
        			itemp4=(int)temp4;
            			phimap3[itemp1-1]=phimap1[itemp2-1];
            			phimap3[itemp3-1]=phimap1[itemp4-1];
			}
         	}
 
		temp=0.0f;
		for( ix=1;ix<=(icgrid-1)/2;ix++) {
			temp=temp + (phimap3[ix-1]*phimap2[ix-1]);
			//System.out.println("phimap3 "+ix+" "+phimap3[ix-1]);
		}
		spec=(2.0f*temp);
		// following needed as spec exceeds 1.0 occasionally in focussing calculations
		// (SS May 8, 1998)
		if(spec>1.0)spec=0.995f;
		System.out.println(" ");
		System.out.println("gauss-seidel spectral radius is"+spec);
		System.out.println(" ");
		return spec;
	}
 

	     void  itit(float spec) {

		//
		//
		// NOTE THIS HAS BEEN ALTERED TO GIVE THE CONVERGENCE TO THE EXACT
		// SOLUTION, WHICH IS INPUTTED VIA OPTION 3,(AND NOW WITH GAUSS
		// SIEDEL IMPLEMENTED, WITH MAPPED SUBARRAYS.)
		//
		//  finally what we`ve all been waiting for-
		// do the actual pb equation iteration
		// first the linear, with relaxation parameter 1.
		// { non-linear with relaxation parameter 0.6
		//
		// this is a modified iteration routine which runs at about
		// three times that of the original version. this has been 
		// accomplished by several features which will be mentioned
		// at the appropiate juncture.
		// note that the arrays slab,old and denom no longer exist
		// and that epsmap is only used in setting up arrays
		// in their place are several arrays, the two biggest being
		// phimap2 and sf.
		// however, the array space accessed in the main inner loop
		// is no larger than before, i.e. memory requirements should
		// not be much different.
		//
		// some notes on array sizes. there can be no more than 10000
		// charges, or 12000 boundary elements, or 40000 high potential salt
		// sites (for the nonlinear). also, for the latter the potential
		// in salt should not exceed 10kt
		//
			 
		int nxran = 60 ;
		int nyran = 20 ;
		float rmsl[]=new float[nxran],rmsn[]=new float[nxran];
		float rmaxl[]=new float[nxran],rmaxn[]=new float[nxran];
	    	float grdn[]=new float[1000];
		String day,symb,title;
		char   iplot[]=new char[nyran];
		boolean qstopper,resdat,once,igt;
		int star,fin,sta1[]=new int[ngrid],sta2[]=new int[ngrid];
		int fi1[]=new int[ngrid],fi2[]=new int[ngrid];
 		int idif1z,iadd2,n,iy,ix,iz,inc1ya,inc2yb,lat2,long1,lat1,itemp2,i;
		float  sixth,temp2,temp1,temp3, temp4,grden;		
		int ihgd2,ihgd,icgrid,iw,je,k,long2,itemp1,itemp4,itemp3,idif1y;
		int inc2xb,inc2xa,idif2x,inc2ya,npoint,ires,inc2za,inc2zb,idif2y,ibin,j,m,midg;
		float rmxch2,rmsch,rnorm2,epsrat,om3,rmxch,res2,res1,ap4,ap3,ap2,ap1,rmsch2,temp;		
		int icount2bbac,ii,nn,itemp,idif2z,inc1yb,inc1za,icount2abac,inc1zb,iadd1,inc1xb,idif1x;
		float res3,res4;
		int inc1xa,itn;
		boolean goto1689,goto1000,goto333;


		ii=0;rmxch2=0;inc2xb=0;grden=0;idif2x=0;inc2ya=0;inc2yb=0;inc2za=0;rmsch2=0;
		idif2y=0;inc2xa=0;idif2z=0;inc2zb=0;inc1xb=0;inc1xa=0;idif1x=0;idif1y=0;inc1yb=0;
		inc1zb=0;inc1za=0;om2=0;idif1z=0;inc1ya=0;
     	 	finish = cputime(finish);
      		day=getTimeStr();
      		System.out.println("now iterating at: "+day);
	
		// allocate memory to arrays

		// some initialization
 
		once=true;
		icount2abac=icount2a;
		icount2bbac=icount2b;
		epsrat=epsout/epsin;
		sixth = 1.f/6.f;
 
		icgrid=igrid*igrid*igrid;
		ihgd=(igrid+1)/2;
		if(icon2==0) {
			icon1=10;
			icon2=1;
		} 
		if(icon1>nlit) icon1=nlit;
 
		for( i = 1;i<=nxran;i++) {
	  		rmsl[i-1] = 0.0f;
	  		rmsn[i-1] = 0.0f;
	  		rmaxl[i-1] =0.0f;
	  		rmaxn[i-1] =0.0f;
		}
        	npoint = (igrid-2)*(igrid-2)*(igrid-2);
 
		// comment out for cray version 14, nov 88, kas
 
		// ---------------------------------------------
		// MAIN SET UP ROUTINE	
		// ---------------------------------------------
		 
		if(iper[0]) {
			n=0;
			for( iz=2;iz<=igrid-1;iz++) {
	  			iadd1=(iz-1)*igrid*igrid ;
	  			for( iy=2;iy<=igrid-1;iy++) {
	    				iadd2=(iadd1+(iy-1)*igrid +2)/2;
	    				n=n+1;
	    				bndx[n-1]=iadd2;
				}
			}
			idif1x=(igrid-2)/2;
			idif2x=idif1x+1;
			inc1xa=1;
			inc1xb=0;
			inc2xa=0;
			inc2xb=1;
		} 
		if(iper[1]) {
			n=0;
			for( iz=2;iz<=igrid-1;iz++) {
	  			iadd1=(iz-1)*igrid*igrid ;
	  			for( ix=2;ix<=igrid-1;ix++) {
	    				iadd2=(iadd1+ix+1)/2;
	   				n=n+1;
    	    				bndy[n-1]=iadd2;
				}
			}
			idif1y=igrid*(igrid-2)/2;
			idif2y=idif1y+1;
			inc1ya=(igrid/2)+1;
			inc1yb=inc1ya-1;
			inc2ya=inc1yb;
			inc2yb=inc1ya;
		}
		if(iper[2]) {
			n=0;
			for( ix=2;ix<=igrid-1;ix++) {
	  			iadd1=ix+1;
	  			for( iy=2;iy<=igrid-1;iy++) {
	    				iadd2=(iadd1+(iy-1)*igrid)/2;
	    				n=n+1;
	    				bndz[n-1]=iadd2;
				}
			}
			idif1z=igrid*igrid*(igrid-2)/2;
			idif2z=idif1z+1;
			inc1za=((igrid*igrid)/2)+1;
			inc1zb=inc1za;
			inc2za=inc1zb;
			inc2zb=inc1za;
		}
 
		// END OF SET UP	
 
		// remove qstopper file if it already exists
 
		/*      inquire[file='qstop.test',exist=qstopper];
      		if(qstopper) {
      			open[30,file='qstop.test');
      			close[30,status='delete');
      			qstopper=false;
      		} 
		*/
 
		// check for resolution data
		// commented out, by kas, 29my 89, as syntax not vax compatible, and obseleted
 
		//      inquire(file='res.dat',exist=resdat)
 
		resdat = false;
      		res1=0.0f;
      		res2=0.0f;
		res3=0.0f;
		res4=0.0f;
 
		if(resdat) { 
			System.out.println(" ");
			System.out.println("linear resolution criteria are:"+res1+" "+res2);
			if(nnit!=0) {
				System.out.println("non-linear resolution criteria are:"+res3+" "+res4);
			}
			System.out.println(" ");
		} 
 
		System.out.println(" ");
		System.out.println(" ");
		if(gten>0.0){
			System.out.println("rms-change     max change    grid energy    #iterations");
		}else {
     			System.out.println("  rms-change     max change       #iterations");
		}
 
		// set up start and stop vectors
		sta1[2-1]=(igrid*igrid + igrid +4)/2;
		sta2[2-1]=sta1[2-1]-1;
		fi1[2-1]=igrid*igrid - (igrid+1)/2;
		fi2[2-1]=fi1[2-1];
		itemp1=igrid + 2;
		itemp2=igrid*igrid -igrid -2;
		for( i=3;i<=igrid-1;i++) {
			sta1[i-1]=fi1[i-1-1] + itemp1;
			sta2[i-1]=fi2[i-1-1] + itemp1;
			fi1[i-1]=sta1[i-1-1] + itemp2;
			fi2[i-1]=sta2[i-1-1] + itemp2;
		}
 
		// also
 
		lat1= (igrid-1)/2;
		lat2= (igrid+1)/2;
		long1= (igrid*igrid - 1)/2;
		long2= (igrid*igrid + 1)/2;
 
      		ires=0;
 		goto333=false;
		if(icheb) {
			om2=1.0f;
			//goto 333;
			goto333=true;
		}
 		if(!goto333) {
		om2=2.0f/(1.0f + (float)Math.sqrt(1 - spec));
		for( ix=1;ix<=(icgrid+1)/2;ix++) {
	  		sf1[ix-1]=sf1[ix-1]*om2;
	  		sf2[ix-1]=sf2[ix-1]*om2;
		}
		for( ix=1;ix<=icount1b ;ix++) {
	  		qval[ix-1]=qval[ix-1]*om2;
		} 
		for( iy=1;iy<=6;iy++) {
		for( ix=1;ix<=icount2b;ix++) {
	  		je=newIndexTwo(ix,iy,6);
	  		db[je]=db[je]*om2;
		}
		}
		sixth=sixth*om2 ;
 		}
		om1=1.0f-om2;
 
		i=1;
 
        	iw=1;
		for( iz=1;iz<=igrid;iz++) {
	 	for( iy=1;iy<=igrid;iy++) {
	   	for( ix=1;ix<=igrid;ix++) {
			je=newIndexThree(ix,iy,iz,igrid);
	    		phimap3[iw-1]=phimap[je];
	    		iw=iw+1;
		}
		}
		}
 		//1689 
		while(true) {
		for( ix=1;ix<=(icgrid+1)/2;ix++) {
			iy=ix*2;
			phimap1[ix-1]=phimap3[iy-1-1];
			phimap2[ix-1]=phimap3[iy-1];
		}
 
		//       set boundary values 
 
		star=(igrid+1)/2;
		iy=(igrid*(igrid+1)/2) - igrid + 1;
		fin=(igrid*(igrid-1)-2)/2;
		ihgd2=ihgd-1;
		for( ix=star;ix<=fin;ix++) {
	  		iy=iy+igrid;
	  		bndx1[ix-1]=phimap1[iy-1];
	  		bndx2[ix-1]=phimap1[iy+ihgd2-1];
		}
  
		star=(igrid+2)/2;
		iy=(igrid*(igrid+2)/2) - igrid +1;
		fin=(igrid*(igrid-1)-1)/2;
		for( ix=star;ix<=fin;ix++) {
     	  		iy=iy+igrid;
	  		bndx3[ix-1]=phimap2[iy-1];
	  		bndx4[ix-1]=phimap2[iy+ihgd2-1];
		}
 
		while(true) {
		//1000      continue
 
		// clear rms, max change
 
	 	rmsch = 0.0f;
	  	rmxch = 0.00f;

 
		// if there is no salt { the main loop is executed without sf
		// saving about 15% in execution time
 
		if(rionst>0.0) {
        		for( n = 2; n<=igrid-1;n++) {
				star=sta1[n-1];
				fin=fi1[n-1];
            			for( ix = star;ix<=fin;ix++) {
            				temp1 = phimap2[ix-1] + phimap2[ix-1-1];
            				temp2 = phimap2[ix+lat1-1] + phimap2[ix-lat2-1];
            				temp3 = phimap2[ix+long1-1] + phimap2[ix-long2-1];
       					phimap1[ix-1] =phimap1[ix-1]*om1 + (temp1+temp2+temp3)*sf1[ix-1];
				}
			}
 
			// otherwise the main loop is as below:
 
        	}else {
 
        		for( n = 2; n<=igrid-1;n++) {
	    			star=sta1[n-1];
	    			fin=fi1[n-1];
            			for( ix = star;ix<=fin;ix++) {
            				temp1 = phimap2[ix-1] + phimap2[ix-1-1];
            				temp2 = phimap2[ix+lat1-1] + phimap2[ix-lat2-1];
            				temp3 = phimap2[ix+long1-1] + phimap2[ix-long2-1];

      					phimap1[ix-1] = phimap1[ix-1]*om1 + (temp1+temp2+temp3)*sixth;
				}
			}
        	}
 
        	// the above loops are about fourtimes faster than the original
        	// loop over all grid points for several reasons, the biggest being that
        	// we are only solving laplace's equation (unless salt is present), which
        	// numerically much simpler, hence faster. we put all we leave out, back
        	// in below, ending up with an equivalent calculation, but much faster.
        	//
        	// first we add back the dielectric boundary points, by recalculating them
        	// individually. note this is still vectorised by means of a gathering
        	// load by the compiler.
        	//
        	//$DIR NO_RECURRENCE 
        	for( k=1;k<=icount2a;k++) {
			ix=idpos[k-1];
			je=newIndexTwo(k,1,6);
			temp1=phimap2[ix-1-1]*db[je]+phimap2[ix-1]*db[je+1];
			temp2=phimap2[ix-lat2-1]*db[je+2]+phimap2[ix+lat1-1]*db[je+3];
        		temp3=phimap2[ix-long2-1]*db[je+4]+phimap2[ix+long1-1]*db[je+5];
        		phimap1[ix-1]= phimap1[ix-1] + temp1+temp2+temp3;
		}
 
		// next we add back an adjustment to all the charged grid points due to
		// the charge assigned. the compiler directive just reassures the vector
		// compiler that all is well as far as recurrence is concerned, i.e. it
		// would think there is a recurrence below, where as in fact there is none.
 
		// Now reset boundary values altered in above loops.
 
		star=(igrid+1)/2;
		iy=(igrid*(igrid+1)/2) - igrid +1;
		fin=(igrid*(igrid-1)-2)/2;
		//$DIR NO_RECURRENCE
		for( ix=star;ix<=fin;ix++) {
	  		iy=iy+igrid;
	  		phimap1[iy-1]=bndx1[ix-1];
	  		phimap1[iy+ihgd2-1]=bndx2[ix-1];
		}
 
		//$DIR NO_RECURRENCE 
		for( k=1;k<=icount1a;k++) {
			temp=qval[k-1];
			ix=iqpos[k-1];
			phimap1[ix-1]=phimap1[ix-1] + temp;
		}
 
		// if periodic boundary condition option
		// force periodicity using wrap around update of boundary values:
		// 2nd slice-->last
		// last-1 slice-->first
 
		// z periodicity
 
        	if(iper[2]) {
          		for( iz = 1;iz<=(igrid-1)*(igrid-1);iz+=2) {
	  			temp1=bndz[iz-1];
	  			temp2=temp1+idif1z;
	  			temp3=temp2+inc1za;
	  			temp4=temp1+inc1zb;
				itemp1=(int)temp1;
				itemp2=(int)temp2;
				itemp3=(int)temp3;
				itemp4=(int)temp4;
	   			phimap1[itemp1-1]=phimap2[itemp2-1];
	    			phimap1[itemp3-1]=phimap2[itemp4-1];
			}
        	}
 
		// y periodicity
 
        	if(iper[1]) {
          		for( iy = 1;iy<=(igrid-1)*(igrid-1);iy+=2) {
	  			temp1=bndy[iy-1];
	  			temp2=temp1+idif1y;
	  			temp3=temp2+inc1ya;
	  			temp4=temp1+inc1yb;
        			itemp1=(int)temp1;
        			itemp2=(int)temp2;
        			itemp3=(int)temp3;
        			itemp4=(int)temp4;
            			phimap1[itemp1-1]=phimap2[itemp2-1];
            			phimap1[itemp3-1]=phimap2[itemp4-1];
			}
        	}
 
		// x periodicity
 
        	if(iper[0]) {
          		for( ix = 1;ix<=(igrid-1)*(igrid-1);ix+=2) {
	  			temp1=bndx[ix-1];
	  			temp2=temp1+idif1x;
	  			temp3=temp2+inc1xa;
	  			temp4=temp1+inc1xb;
        			itemp1=(int)temp1;
        			itemp2=(int)temp2;
        			itemp3=(int)temp3;
        			itemp4=(int)temp4;
           	 		phimap1[itemp1-1]=phimap2[itemp2-1];
            			phimap1[itemp3-1]=phimap2[itemp4-1];
			}
         	}
 
		if(icheb) {
       		 	//**omalt(sf1,sf2,qval,db,sixth,2*i-1,om1,om2,spec,rionst ,nhgp,icount1b,icount2b];
			//omalt(sixth,2*i-1,spec);
			sixth=omalt(sixth,2*i,spec,rionst,nhgp,icount1b,icount2b);
		}
 
		// Next update phimap2 using the new phimap1
 
        	if(rionst>0.0) {	
        		for( n = 2; n<=igrid-1;n++) {
	    			star=sta2[n-1];
	    			fin=fi2[n-1];
            			for( ix = star;ix<=fin;ix++) {
            				temp1 = phimap1[ix-1] + phimap1[ix+1-1];
 
            				temp2 = phimap1[ix+lat2-1] + phimap1[ix-lat1-1];
 
            				temp3 = phimap1[ix+long2-1] + phimap1[ix-long1-1];
 
       	    				phimap2[ix-1] =phimap2[ix-1]*om1 +  (temp1+temp2+temp3)*sf2[ix-1];
				}
			}
		}
 
		else {
 
        		for( n = 2; n<=igrid-1;n++) {
	    			star=sta2[n-1];
	    			fin=fi2[n-1];
            			for( ix = star;ix<=fin;ix++) {
            				temp1 = phimap1[ix-1] + phimap1[ix+1-1];
 
            				temp2 = phimap1[ix+lat2-1] + phimap1[ix-lat1-1];
 
            				temp3 = phimap1[ix+long2-1] + phimap1[ix-long1-1];
 
       	    				phimap2[ix-1] =phimap2[ix-1]*om1 + (temp1+temp2+temp3)*sixth;
				}
			}
		}
 
		//$DIR NO_RECURRENCE 
        	for(k=icount2a+1;k<=icount2b;k++) {
			ix=idpos[k-1];
			je=newIndexTwo(k,1,6);
			temp1=phimap1[ix-1]*db[je]+phimap1[ix+1-1]*db[je+1];
			temp2=phimap1[ix-lat1-1]*db[je+2]+phimap1[ix+lat2-1]*db[je+3];
        		temp3=phimap1[ix-long1-1]*db[je+4]+phimap1[ix+long2-1]*db[je+5];
        		phimap2[ix-1]=phimap2[ix-1] + temp1+temp2+temp3;
		}
		// reset boundary condition
 
		star=(igrid+2)/2;
		iy=(igrid*(igrid+2)/2) - igrid +1;
		fin=(igrid*(igrid-1)-1)/2;
		//$DIR NO_RECURRENCE
		for( ix=star;ix<=fin;ix++) {
			iy=iy+igrid;
			phimap2[iy-1]=bndx3[ix-1];
			phimap2[iy+ihgd2-1]=bndx4[ix-1];
		}
 
		//$DIR NO_RECURRENCE 
		for( k=icount1a+1;k<=icount1b;k++) {
			temp=qval[k-1];
			ix=iqpos[k-1];
			phimap2[ix-1]=phimap2[ix-1] + temp;
		}
 
		// z periodicity
 
        	if(iper[2]) {
          		for( iz = 2;iz<=(igrid-1)*(igrid-1);iz+=2) {
	  			temp1=bndz[iz-1];
	  			temp2=temp1+idif2z;
	  			temp3=temp2+inc2za;
	  			temp4=temp1+inc2zb;
        			itemp1=(int)temp1;
        			itemp2=(int)temp2;
        			itemp3=(int)temp3;
        			itemp4=(int)temp4;
            			phimap2[itemp1-1]=phimap1[itemp2-1];
            			phimap2[itemp3-1]=phimap1[itemp4-1];
			}
       		}
 
		// y periodicity
 
        	if(iper[1]) {
           		for(iy = 2;iy<=(igrid-1)*(igrid-1);iy+=2) {
	  			temp1=bndy[iy-1];
	  			temp2=temp1+idif2y;
	  			temp3=temp2+inc2ya;
	  			temp4=temp1+inc2yb;
        			itemp1=(int)temp1;
        			itemp2=(int)temp2;
        			itemp3=(int)temp3;
        			itemp4=(int)temp4;
            			phimap2[itemp1-1]=phimap1[itemp2-1];
            			phimap2[itemp3-1]=phimap1[itemp4-1];
			}
        	}
 
		// x periodicity
 
        	if(iper[0]) {
          		for( ix = 2;ix<=(igrid-1)*(igrid-1);iy+=2) {
	  			temp1=bndx[ix-1];
	  			temp2=temp1+idif2x;
	  			temp3=temp2+inc2xa;
	  			temp4=temp1+inc2xb;
        			itemp1=(int)temp1;
        			itemp2=(int)temp2;
        			itemp3=(int)temp3;
        			itemp4=(int)temp4;
           		 	phimap2[itemp1-1]=phimap1[itemp2-1];
            			phimap2[itemp3-1]=phimap1[itemp4-1];
			}
         	}
 
		if(icheb) {
        		//**omalt(sf1,sf2,qval,db,sixth,2*i,om1,om2,spec,rionst ,nhgp,icount1b,icount2b];
			sixth=omalt(sixth,2*i,spec,rionst,nhgp,icount1b,icount2b);
			
		}
 
		// we also save time by only checking convergence every ten
		// iterations, rather than every single iteration.
 
		if(mod(i,icon1)==(icon1-1)) {
			for( ix=2;ix<=(icgrid+1)/2;ix+=icon2) {
	  			phimap3[ix-1]=phimap2[ix-1];
			}
		}
 
		if(gten>0.0) {
			grden=0.0f;
			for(ix=1;ix<=icount1a;ix++) {
				iy=iqpos[ix-1];
				grden=grden + phimap1[iy-1]*gval[ix-1];
			}
			for(ix=icount1a+1;ix<=icount1b;ix++) {
				iy=iqpos[ix-1];
				grden=grden+phimap2[iy-1]*gval[ix-1];
			}
			grdn[i-1]=grden/2.0f;
			if(i>10) {
				igt=true;
				for(ix=i-4;ix<=i;ix++) {
				for( iy=i-4;iy<=i;iy++) {
					if(Math.abs(grdn[iy-1]-grdn[ix-1])>gten) igt=false;
				}
				}
				if(igt){
					for(iy=i-4;iy<=i;iy++) {
						System.out.println(grdn[iy-1]);
					}
					ires=1;
				}
			}
 
		}
 
		if((mod(i,icon1)==0)||(ires==1)) {
			rnorm2=0;
        		for( ix=2;ix<=(icgrid+1)/2;ix+=icon2) {
	    			temp2=phimap3[ix-1]-phimap2[ix-1];
	    			rnorm2=rnorm2+temp2*temp2;
	    			rmxch=Math.max(rmxch,Math.abs(temp2));
			}
        		rmsch = (float)Math.sqrt((float)(icon2)*rnorm2/npoint);
			rmsch2=rmsch;
			rmxch2=rmxch ;
			if(gten>0.){
	 			System.out.println(rmsch2+" "+rmxch2+" "+grden + "at "+i+"iterations");
			}else {
				System.out.println(  rmsch2+ " "+rmxch2+" at"+i+"iterations");
			}
          		if(rmsch<=res1&&rmxch<=res2) ires=1;
	    		if(igraph&&(once)) { 
	    			for( j=i-9;j<=i;j++) {
	    				ibin =(int)( (j-1.f)*(nxran-1.f)/(nlit-1.f) + 1);
	    				rmsl[ibin-1] = rmsch;
           				rmaxl[ibin-1] = rmxch;
				}
	    		}
        	}
 
		// check to see if accuracy is sufficient, or if a qstop command
		// has been issued
 
      
 
      		i=i+1;
 		goto1000=false;
		if(gten<1.e-6) {
			if(i<=nlit&&ires==0) goto1000=true; //goto 1000;
		}else {
			if(ires==0) goto1000=true; //goto 1000;
		}
		if(goto1000) continue;
		else 	     break;
 		}
		//	end of iteration loop
 
		//       remap into phimap
 
		for( iy=1;iy<=(icgrid-1)/2 ;iy++) {
			ix=iy*2 ;
			phimap3[ix-1-1]=phimap1[iy-1];
			phimap3[ix-1]=phimap2[iy-1];
		}
		if(once) {
			iw=1;
			for( iz=1;iz<=igrid;iz++) {
	  		for( iy=1;iy<=igrid;iy++) {
	    		for( ix=1;ix<=igrid;ix++){
				je=newIndexThree(ix,iy,iz,igrid);
	       			phimap[je]=phimap3[iw-1];
	     			iw=iw+1;
			}
			}
			}
			je=newIndexThree(igrid,igrid,igrid,igrid);
			phimap[je]=phimap1[(icgrid+1)/2-1];
		}else {
			iw=1 ;
			for( iz=1;iz<=igrid;iz++){
	 		for( iy=1;iy<=igrid;iy++) {
	   		for( ix=1;ix<=igrid;ix++) {
				je=newIndexThree(ix,iy,iz,igrid);
	       			phimap[je]=phimap[je]-phimap3[iw-1];
	     			iw=iw+1;
			}
			}
			}
			je=newIndexThree(igrid,igrid,igrid,igrid);
			phimap[je]=phimap[je]  - phimap1[(icgrid+1)/2-1];
		} 
 		goto1689=false;
		if(diff&&(once)) {
			System.out.println("now doing uniform dielectric run....");
			once=false;
			om3= 2.0f/(1 +(3.14159f/(float)(igrid)));
			sixth=sixth*om3/om2;
			for( i=1;i<=icgrid+1;i++) {
				phimap3[i-1]=phimap3[i-1]*epsrat ;
			}
			for( i=1;i<=icount1b;i++) {
				qval[i-1]=qval[i-1]*om3/om2;
			}
			for( i=1;i<=ibc;i++) {
				je=newIndexTwo(i,5,5);
				itemp=(int)cgbp[je];
				je=newIndexTwo(i,4,5);
				qval[itemp-1]=cgbp[je]*om3;
			}
			om1=1-om3 ;
			icount2a=0;
			icount2b=0;
			nlit=(int)(7.8f*(float)(igrid)/3.14159f);
			i=1;
			goto1689=true;//goto 1689;
		}
		if( goto1689) continue;
		else break;
 		}
		icount2a=icount2abac;
		icount2b=icount2bbac;
 
      		finish = cputime(finish);
        	day=getTimeStr();
		System.out.println("finished qdiffx linear iterations");
		System.out.println("at                       : "+day);
        	System.out.println("total time elapsed so far: "+finish);
		System.out.println("# loops                  : "+(i-1));
		System.out.println("mean,max change (kT/e]   :"+rmsch2+" "+rmxch2);
 
		//  plot convergence history   
  
		// code phimap corner, for use in transference from irises to convex
		// and via versa
		je=newIndexThree(1,1,1,igrid);
		ap1=phimap[je];
		ap2=ap1*10000;
		ap3=(float)((int)(ap2));
		if(ap3>0) {
			ap4=(ap3+0.8f)/10000;
		}else {
			ap4=(ap3-0.8f)/10000;
		}
	
		phimap[je]=ap4;
 
		if(ipoten) {
			midg = (igrid+1)/2;
			for( m = 1;m<=5;m++) {
	  			n = (igrid - 1)/4;
	  			nn = (m-1)*n + 1;
	  			System.out.println("phi"+nn+" "+midg);
				je=newIndexThree(nn,midg,ii,igrid);	
	  			for(ii=0;ii<=igrid;ii++) {
					je=newIndexThree(nn,midg,ii,igrid);	
	  				System.out.println((phimap[je])); 
				}
			}
		} 

	}

 

	
	float  omalt(float sixth,int it,float spec,float salt,int nhgp,int ncrg,int nsp) {
  	
	        float om4,om3;
		int je,iy,ix,iz;
		om3=1/(1-om2*spec*0.25f);
 		 
		 
	
		if(om1<1.e-6) om3=1.f/(1.f-om2*spec*0.5f);
 
		om4=om3/om2;
 
		om2=om3;
		om1=1.0f-om2;
 
		if(salt>0.0f) {
			if(mod(it,2)==1) {
        			for( ix=1;ix<=nhgp;ix++) {
          				sf1[ix-1]=sf1[ix-1]*om4;
				}
			}else {
				for( ix=1;ix<=nhgp;ix++) {
          				sf2[ix-1]=sf2[ix-1]*om4;
				}
			}
		}
 
        	for( ix=1;ix<=ncrg;ix++) {
          		qval[ix-1]=qval[ix-1]*om4;
		}
 
        	for( iy=1;iy<=6;iy++) {
        	for( ix=1;ix<=nsp;ix++) {
			je=newIndexTwo(ix,iy,6);
          		db[je]=db[je]*om4;
		}
		}
 
        	sixth=sixth*om4;
  		return sixth;
	}
 
	public static void main(String argv[]) {

		StructFileReader ss=new StructFileReader("2.pdb");
		ss.setAtomCoo();
		Molecule molecule = ss.getMolecule();
		for(Molecule m=molecule;m!=null;m=m.next)
                for(Strand s=m.strand;s!=null;s=s.next)
                for(Residue r=s.residue;r!=null;r=r.next)
                for(Atom a=r.atom;a!=null;a=a.next) {
                      	a.radius=2.f;
			//if(a.name.startWith("N")) 
			a.charge=1.f;
                }
		
		Delphi delphi=new Delphi(molecule);
		delphi.logs=true;
		delphi.isrf=true;
		delphi.setPrm(65,70.0f, 2.f,80.f,1.6f,0f,0.14f);
		
		delphi.qdiff();
	}

	

	void setcrg(){

		//
		// gchrg is the fractional charge in electron units assigned to
		// each grid point, gchrgp is the position of each such charge on the
		// grid. 
						
		boolean isum[]=new boolean[3*ngrid];
		float gchrgd[];		
		int gchrg2[]=new int[ngcrg];
		int ix,iy,iz,je,ke,kb2,ig,j,kb1,i,k,kb3,n,ico,i2,ib,itemp,ieps,i1,iv,iw,isgrid;
		float debfct,fpoh,sixeps,difeps,cg3,cg2,cg1,epsins6,epsval;

		// assign fractional charges to grid points, using phimap
		// as a dummy array (set to zero again later)

        	for( i=1;i<=igrid;i++) {
        	for( j=1;j<=igrid;j++) {
        	for( k=1;k<=igrid;k++) {
			je=newIndexThree(k,j,i,igrid);
        	 	phimap[je]=0.0f;
       	 	}
        	}
        	}
  
		difeps=epsin-epsout;
		sixeps=epsout*6.0f;
		fpoh=fpi*scale;
		debfct = epsout/((deblen*scale)*(deblen*scale));
 
		for( ig=1;ig<=nqgrd;ig++) {
			je=newIndexTwo(ig,1,4);
			kb1=(int)chrgv2[je];
			kb2=(int)chrgv2[je+1];
			kb3=(int)chrgv2[je+2];
			//System.out.println("nqgrd "+ig+" "+chrgv2[je]+" "+chrgv2[je+1]+" "+chrgv2[je+2]);
	    		for( ix=0;ix<=1;ix++) {
	    			i=kb1+ix;
				je=newIndexTwo(ig,1,4);
	    			cg1=kb1-chrgv2[je]+1-ix;
	      			for(iy=0;iy<=1;iy++){
	      				j=kb2+iy;
					je=newIndexTwo(ig,2,4);
	      				cg2=kb2-chrgv2[je]+1-iy;
					for( iz=0;iz<=1;iz++) {
						k=kb3+iz;
						je=newIndexTwo(ig,3,4);
						cg3=kb3-chrgv2[je]+1-iz;
						ke=newIndexThree(i,j,k,igrid);
						je=newIndexTwo(ig,4,4); 
						//System.out.println("phimap "+ke+" "+je+" "+ig+" "+i+" "+j+" "+k);
	    					phimap[ke]=phimap[ke]+Math.abs(cg1*cg2*cg3)*chrgv2[je];
					}
				}
			}
		}
 

		n=0;
        	for(k=2;k<=igrid-1;k++) {
          		for(j=2;j<=igrid-1;j++){
            			for(i=2;i<=igrid-1;i++){
					je=newIndexThree(i,j,k,igrid);
                			if(phimap[je]!=0) n=n+1;
	    			}
	  		}
		}
 
		gchrgp=new int[3*n];
		gchrg=new float[n];
		gchrgd=new float[n];
		gchrg2=new int[n];
		qval=new float[n];
		iqpos=new int[n];
		gval=new float[n];
		// set up odd/even boolean array 
		for( i=1;i<=3*igrid;i+=2) isum[i-1]=true;
		for( i=2;i<=(3*igrid-1);i+=2)  isum[i-1]=false;
 
		// find which grid points have charge assigned to them
		// (will use this array later to calculate grid energy)
  
		n=0;
		for( k=2;k<=igrid-1;k++) {
	  	for( j=2;j<=igrid-1;j++) {
	    	for( i=2;i<=igrid-1;i++) {
			ke=newIndexThree(i,j,k,igrid);
			if(phimap[ke]!=0) {
				n=n+1;
				je=newIndexTwo(n,1,3);
				gchrgp[je]=i;
				gchrgp[je+1]=j;
				gchrgp[je+2]=k;
				gchrg[n-1]=phimap[ke];
				phimap[ke]=0.0f;
			}
		}
		}
		}
		icount1b=n;
 
		if(iwgcrg) {
			epsval=repsin;
			//call wrtgcrg(gcrgnam,gcrglen,icount1b,gchrgp,gchrg,epsval,scale);
 
		}
		// determine how many charged grid points are odd
 
		icount1a=0;
		for( i=1;i<=n;i++) {
			je=newIndexTwo(i,1,3);
			itemp=gchrgp[je]+gchrgp[je+1]+gchrgp[je+2];
			if(isum[itemp-1]) icount1a=icount1a+1;
		}
 
		// set up odd/even pointer array, to be used in making qval
		// and iqpos
  
		i1=0;
		i2=icount1a;
		for( i=1;i<=n;i++) {
			je=newIndexTwo(i,1,3);
			itemp=gchrgp[je]+gchrgp[je+1]+gchrgp[je+2];
			if(isum[itemp-1]) {
				i1=i1+1;
				gchrg2[i-1]=i1;
			}else {
				i2=i2+1;
				gchrg2[i-1]=i2;
			}
		}
 
		// determine denominator at all charged grid points
 
		ib=0;
		ico=0;
		epsins6=sixeps+ 6.0f*difeps ;
		
		cgbp=new float[5*n+10];
		for( i=1;i<=n;i++) {
			je=newIndexTwo(i,1,3);
			iz=gchrgp[je+2];
			iy=gchrgp[je+1];
			ix=gchrgp[je];
			ieps=0;
			je=newIndexFour(ix,iy,iz,1,igrid);
			if(iepsmp[je]!=0) ieps=ieps+1;
			if(iepsmp[je+1]!=0) ieps=ieps+1;
			if(iepsmp[je+2]!=0) ieps=ieps+1;
			je=newIndexFour(ix-1,iy,iz,1,igrid);
			if(iepsmp[je]!=0) ieps=ieps+1;
			je=newIndexFour(ix,iy-1,iz,2,igrid);
			if(iepsmp[je]!=0) ieps=ieps+1;
			je=newIndexFour(ix,iy,iz-1,3,igrid);
			if(iepsmp[je]!=0) ieps=ieps+1;
			itemp=ieps;
			je=newIndexThree(ix,iy,iz,igrid);
			gchrgd[i-1]=itemp*difeps + debfct*idebmap[je] + sixeps ;
			if(itemp!=6) {
				ib=ib+1;
				if(itemp==0) {
					ico=ico+1;
				}
				je=newIndexTwo(ib,1,5);
				ke=newIndexThree(ix,iy,iz,igrid);
				//System.out.println("cgbp...."+i+" "+ib+" "+ix+" "+iy+" "+iz+" "+n);  
				cgbp[je]=  ix;
				cgbp[je+1]=iy;
				cgbp[je+2]=iz;
				cgbp[je+3]=gchrg[i-1]*fpoh/(epsins6+debfct*idebmap[je]);
				cgbp[je+4]=gchrg2[i-1];
			} 
			//System.out.println("gchrgd "+i+" "+gchrgd[i-1]+" "+ix+" "+iy+" "+iz);
		}
		ibc=ib;
		System.out.println("no. grid points charged and at boundary="+ibc);
		if(ico!=0) System.out.println("## "+ico+"charges are in solution ##");
 
		// make qval, fpoh term so potentials will be in kt/e
 
		for( i=1;i<=n;i++) {
			j=gchrg2[i-1];
			qval[j-1]=gchrg[i-1]*fpoh/gchrgd[i-1];
			gval[j-1]=gchrg[i-1];
		}
 
		// make iqpos
 
		isgrid=igrid*igrid;
		for( i=1;i<=n;i++) {
			j=gchrg2[i-1];
			je=newIndexTwo(i,1,3);
			ix=gchrgp[je];
			iy=gchrgp[je+1];
			iz=gchrgp[je+2];
			iw=1+ix+igrid*(iy-1)+isgrid*(iz-1);
			iv=iw/2;
			iqpos[j-1]=iv;
		}
 
		// end of chrgup, return with qval,iqpos and gchrgp and gchrg
		// also icount1a, icount1b 
		
	}
 
	
	


        float phintp(float gp[]){
		//
		// interpolates the potential at any point inside
		// a cubical volume using the potential values at the
		// 8 vertices by means of a trilinear function:
		//	
		// W = A1.XYZ + A2.XY + A3.XZ + A4.YZ + A5.X +
		//		A6.Y + A7.Z + A8
		//	
		// where Ai coefficients are linear combinations of Wi at
		// the cube corners
		//
	 
 		float phi,rgrid,xg,yg,zg,xgr,ygr,zgr,a8,a7,a6,a5,a4,a3,a2,a1;
      		int nx,ny,nz,nx1,ny1,nz1,k;
 
		// return 0.0 if outside grid
 
		rgrid = (float)(igrid);
		for( k = 1;k<=3;k++) {
	  		if((gp[k-1]<1)||(gp[k-1]>rgrid)) {
	    			phi = 0.0f;
	   			return phi;
	  		}
 		}
 
 
		// find lower left bottom grid point
 
      		nx = (int)(gp[0]);
		nx1 = nx+1;
		if(nx1>igrid)nx1=nx;
      		ny = (int)(gp[1]);
		ny1 = ny+1;
		if(ny1>igrid)ny1=ny;
      		nz = (int)(gp[2]);
		nz1 = nz+1;
		if(nz1>igrid)nz1=nz;
 
 		// calculate cube coordinates of point
 
		xg = nx;
        	yg = ny;
		zg = nz;
		xgr = gp[0] - xg;
		ygr = gp[1] - yg;
		zgr = gp[2] - zg;
 
		// calculate coefficients of trilinear function
 		 
		a8 = phimap[newIndexThree(nx,ny,nz,igrid)];
		a7 = phimap[newIndexThree(nx,ny,nz1,igrid)] - a8;
		a6 = phimap[newIndexThree(nx,ny1,nz,igrid)] - a8;
		a5 = phimap[newIndexThree(nx1,ny,nz,igrid)] - a8;
		a4 = phimap[newIndexThree(nx,ny1,nz1,igrid)] - a8 - a7 - a6;
		a3 = phimap[newIndexThree(nx1,ny,nz1,igrid)] - a8 - a7 - a5;
		a2 = phimap[newIndexThree(nx1,ny1,nz,igrid)] - a8 - a6 - a5;
		a1 = phimap[newIndexThree(nx1,ny1,nz1,igrid)] - a8 - a7 - a6 - a5 - a4 - a3 - a2;
 
		// determine value of phi
		 
		phi = a1*xgr*ygr*zgr + a2*xgr*ygr + a3*xgr*zgr + a4*ygr*zgr + a5*xgr + a6*ygr + a7*zgr + a8;

     		return phi;
	}

	void gtoc(float g[], float c[]){
 
		float goff = (igrid + 1.f)/2.f;
      		for(int i = 1;i<=3;i++) {
        		 c[i-1] = (g[i-1] - goff)/scale + oldmid[i-1];
       		}
	}

 
      	void ctog(float c[],float g[]){
 
		// converts float to grid coordinates
 
	 
		float goff = (igrid + 1.f)/2.f;
      		for(int i = 1;i<=3;i++) {
         		g[i-1] = (c[i-1] - oldmid[i-1])*scale + goff;
       		}
	}
	

	void encalc(){
 
		float epsval=epsin*epkt; 
		react(epsval); 
	}

	void react(float epsval){
 
 		float temp1,temp2,temp3,sixth,fact;
		int i,je,ix,iy,iz;

		//schrg= new float[ibnum];
		schrg= new float[surface.vtot]; 
		
        	
		// fact=six divide by 4 pi, multiplied by epkt, scale, divided by epsval
 
		fact=epsval*0.9549296586f/(2.0f*scale*epkt);
		sixth=1.0f/6.0f;
		 
 		surface.resetVerticesNumber();
		surface.calcVerticeArea();

		// calculate surface charge
 		float ent=0;
		 
		float coo[]=new float[3];
		float goo[]=new float[3];
		for( i=0;i<surface.vtot;i++){
			
			je=i*3;
			coo[0]=surface.vert[je];
			coo[1]=surface.vert[je+1];
			coo[2]=surface.vert[je+2];
			ctog(coo,goo); 
			ix=(int) (goo[0]+0.5);
			iy=(int) (goo[1]+0.5);
			iz=(int) (goo[2]+0.5);
			 
	 		
        		temp1=phimap[newIndexThree(ix+1,iy,iz,igrid)]+phimap[newIndexThree(ix-1,iy,iz,igrid)];
        		temp2=phimap[newIndexThree(ix,iy+1,iz,igrid)]+phimap[newIndexThree(ix,iy-1,iz,igrid)];
        		temp3=phimap[newIndexThree(ix,iy,iz+1,igrid)]+phimap[newIndexThree(ix,iy,iz-1,igrid)];
        		schrg[i]=phimap[newIndexThree(ix,iy,iz,igrid)]-(temp1+temp2+temp3)*sixth;
			//schrg[i-1]=schrg[i-1]*surface.area[i-1]*fact;
			//schrg[i-1]=schrg[i-1]*fact;
			ent+=schrg[i];
			//System.out.println("schrg "+i+" "+schrg[i]);
        	}
		 
		surface.schrg=schrg;
  		System.out.println("total surface induced charges: "+ent);
		// the following notes are no longer relevant. schrg now
		// contains "float" charges.
		// NB float charge is schrg*epsval/(scale*epkt), epkt =561
		// epsval = inner dielectric
		// convert to 'float' charges with fact, these 'charges {
		// carry the dielectric, and scale factors with them
 
		// attempt to spread charge about a bit..
 
		float dist1,dist2,dist3,dist; 	
		float sdist1,sdist2,sdist3;
		float ptemp;
		int   nscrg=ibnum,j;
		
		spot=  new float[surface.vtot];  
		System.out.println("nscrg "+nscrg+" "+surface.vtot);
		 
		for( i=0;i<surface.vtot;i++) {
			
			je=i*3;
			 
			coo[0]=surface.vert[je];
			coo[1]=surface.vert[je+1];
			coo[2]=surface.vert[je+2];
			 
			ctog(coo,goo);
			//ptemp=phintp(goo)*25.6f/1000;
			ptemp=phintp(goo);
			spot[i]=ptemp;
			//System.out.println("spot  "+i+" "+spot[i]);					
		}
 		surface.spot=spot;
		//surface.calPotentialCol();
	}
}
