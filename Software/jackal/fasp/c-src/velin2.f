
	character*31 line
	character*3 atom
	character*80 check
	open(10,file='acc.res')
	open(11,file='temp.pdb')
50      read(11,'(a80)',end=100) check
c       if(check(1:6).ne.'COMPND') go to 50
c       if((check(11:13).eq.'DNA').or.(check(11:13).eq.'RNA')) go to 100
	read(10,'(a31)') line
	i=0
	i1=0
	i2=0
	a=0
	b=0
20      read(10,'(a31)',end=30) line
	if(line(2:6).eq.'total') go to 30
	i=i+1
	read(unit=line(20:27),'(f7.4)') acc
	read(unit=line(3:5),'(a3)') atom
        if(atom.eq.'PHE') a=a+(0-(acc/325.))*3.7
	if(atom.eq.'MET') a=a+(0-(acc/309.))*3.4
	if(atom.eq.'ILE') a=a+(0-(acc/280.))*3.1
	if(atom.eq.'LEU') a=a+(0-(acc/290.))*2.8
	if(atom.eq.'VAL') a=a+(0-(acc/263.))*2.6
        if(atom.eq.'CYS') a=a+(0-(acc/243.))*2.0
	if(atom.eq.'TRP') a=a+(0-(acc/360.))*1.9
c       if(atom.eq.'ALA') a=a+(0-(acc/218.))*1.6
        if(atom.eq.'THR') a=a+(0-(acc/254.))*1.2
        if(atom.eq.'GLY') a=a+(0-(acc/190.))*1.0
        if(atom.eq.'SER') a=a+(0-(acc/235.))*0.6
c       if(atom.eq.'PRO') a=a+(0-(acc/249.))*-0.2
c       if(atom.eq.'TYR') a=a+(0-(acc/327.))*-0.7
c       if(atom.eq.'HIS') a=a+(0-(acc/301.))*-3.0
c       if(atom.eq.'GLN') a=a+(0-(acc/296.))*-4.1
c       if(atom.eq.'ASN') a=a+(0-(acc/271.))*-4.8
c       if(atom.eq.'GLU') a=a+(0-(acc/301.))*-8.2 
c       if(atom.eq.'LYS') a=a+(0-(acc/310.))*-8.8
c       if(atom.eq.'ASP') a=a+(0-(acc/269.))*-9.2
c       if(atom.eq.'ARG') a=a+(0-(acc/345.))*-12.3
	go to 20
30      continue
c       if(check(73:76).eq.'    ') go to 100
c       if((ah.lt.0.15).or.(ap.lt.0.15)) go to 100
	write(*,'(a4,4x,i6,4x,f8.3)') check(73:76),i,a/i
100     continue
	stop
	end
