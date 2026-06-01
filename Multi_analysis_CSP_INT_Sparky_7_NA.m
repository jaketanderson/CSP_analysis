%This script was written by Dominique Frueh and all rights to it are reserved by him

%CSP analysis multiple points, subplots with both absolute and local
%scaling
%requires consistent Sparky lists for each point
%basic analysis for CSP
%provides bar graph, mean, std, median
%the graph can be used to decide which selection to use for CSP
% e.g. all residues above mean + 1std, mean+2std, median+1std etc.
%residue titration displayed as sub-plots
% has residues in subplots 2 by 3, with both scale adjusted to max shift and autoscale
%this allows to visualy judge how important the shifts are.

%update/tailored for ArCP to include signal intensities
%reads Sparky with format res# N(ppm) H(ppm) Int

%CSP now calculated on the fly in matlab script
%use alpha scaling factor (SW(H)/SW(N), or gyromagnetic, or
%residue-corrected
%most studies always employ 1/5 or 0.1 as default

%intensities analysed as McAllister-Driscoll (A-X)/(A-N)
%However, this does not give the graph they plot (you would get build-up,
%they get decay)
%also do (A-X)/A which gives a decay

%two ways or normalizing spectra for losses: to be commented/uncommented
% (i) provide scaling factors measured through DSS, requires DSS at same
% concentration in both titrant and 15N solute OR (bad) calculate dilution
% of DSS if DSS only in 15N solute.
%(ii) normalize to largest signal in spectrum: assumes a C-term or N-term
%residue that does not interact and hence is always the strongest. If using
%ii, also check the residue# at max is the same for all.

'requires shifts and int for all selected residues at all points, bad residues tagged in code'
'Sparky peak lists: res#, N (ppm) H(ppm) Int; recommended: use file named as shift_HNXX with XX = percentage of titrant'

'FIX in progess: apply residue offsets as needed (solution provided untested), likely loopholes'

clear all
close all

%USER can modify ___________________________
%alpha=1/5;          %nitrogen scaling factor (Williamson PNMRS 2015): 0.1-0.45 (1/5 default)

res_off=0;          %allows for changing numbering between constructs or to match full length
shiftY=0.14;           %scaling Y axis in files overallXX.eps (change this for comparing between different titrations with same Y axis)
shiftslopeY=0;      %scaling Y axis in top panel of file CSP_INT_stacked.eps (change this for comparing between different titrations with same Y axis)
intslopeY=0;        %scaling Y axis in bottom panel of file CSP_INT_stacked.eps(change this for comparing between different titrations with same Y axis)
%purging bad residues
%work around for not recoding, should be fixed
residue_bad=[];
%residue_bad=[3;15;18;33;44;45;48;49;50;70;72;81;87;100;104;123;124;125;130;131;132;133;134;135;153;154;159;196;197;254;255;256;257;258;259;260;273;330;356;363;368;371;390;416;417;420;434]-res_off;
%residue_bad=[35;56;70;71;83;120;131;139]
%end USER can modify ___________________________


sc=[];
scalsp=input('do you have scaling factors for each spectrum (y/n)','s');
n=input('number of points in titration? include 0');
Sample=input('Sample name (will be used as a title; keep short)','s');
nc=ceil(n/3);   %used for sub-plots
ALLSHIFT=[];
SHIFT=[];
SHIFTH=[];
SHIFTN=[];
INT=[];
percentage=[];
Ptot=input('concentration of 15N protein')
resrange=[0 80];       %used for seting x-axis, use when comparing various analyses.

for m=1:n   %build array; size = #residues X #titration points

    PC=input('percentage titrant, one digit (e.g 25% enter 25.0)? **start from 0**');
    filename=input('sparky peak list name (format res# Nppm Hppm Int)','s');
    IN=load(filename, '-ASCII');
    if (m == 1)
        REF=IN
        residue=IN(:,1)-res_off;
    end
    
percentage=[percentage PC];

SHIFTN=[SHIFTN (IN(:,2)-REF(:,2))];
SHIFTH=[SHIFTH (IN(:,3)-REF(:,3))];

INT=[INT IN(:,4)];
if (scalsp(1)=='y')
    sc(m)=input('spectrum scaling factor');
end


end

%SHIFTN=alpha*SHIFTN;
shift_good=[];
shiftH_good=[];
shiftN_good=[];
residue_good=[];
INT_good=[];


 %nul=zeros(0,1);
a=0;
j=0;
 for i=1:size(residue,1)
     %find(residue_bad==residue_purge(i,1))
     if (size(find(residue_bad==residue(i,1)),1) == 0)
         j=j+1;
         residue_good(j,:)=residue(i,1);
         %shift_good(j,:)=SHIFT(i,:);
         shiftH_good(j,:)=SHIFTH(i,:);
         shiftN_good(j,:)=SHIFTN(i,:);
         INT_good(j,:)=INT(i,:);
     else
         a=a+1;
     end
 end
 
 residue_purge=residue_good;
 SHIFTH_purge=shiftH_good;
 SHIFTN_purge=shiftN_good;
 INT_purge=INT_good;
 
 %residue=residue_good;
 %SHIFT=shift_good;
 %SHIFTH=shiftH_good;
 %SHIFTN=shiftN_good;
 %INT=INT_good;
 %end work arounf or not recoding
 

SHIFT_purge=((1/2)*(SHIFTH_purge.^2+0.2*(SHIFTN_purge.^2))).^0.5;



[Mx,Imx]=max(INT_purge);

if (size(sc))
    INTS=INT_purge.*((sc')*ones(1,size(INT_purge,1)))';
    INT_purge=INTS;
elseif ((Imx./max(Imx)) == ones(1,size(INT_purge,2)))
    mxres=input(['residue ' int2str(Imx(1)) ' may be used for scaling;scale? (y/n)'],'s');
    if (mxres(1) == 'y')
    INTS=INT./((Mx')*ones(1,size(INT_purge,1)))';
    INT_purge=INTS;
    end
else
    'no scaling performed, your maxima occur for different residues during titration'
end

INT_A=INT_purge(:,1)*ones(1,size(INT_purge,2));
INT_B=INT_purge(:,size(INT_purge,2))*ones(1,size(INT_purge,2));

D_INT=(INT_A-INT_purge)./(INT_A-INT_B);       %McAllister-Driscoll Biochemistry 1996
N_INT=INT_purge./(INT_A);
DN_INT=(INT_A-INT_purge)./(INT_A);

%SHIFT_purge=SHIFT(find(SHIFT>0));
%residue_purge=residue(find(SHIFT>0));


 
%bbb=input('egfseg')    %for debugging
% SHIFT_med=median(SHIFT_purge);
% SHIFT_mean=mean(SHIFT_purge);
% SHIFT_std=std(SHIFT_purge);
 
 %SHIFT_med=[SHIFT_med SHIFT_med];
 
 %{
     
 figure(m)
 clf
 bar(residue_purge,SHIFT_purge)

 hold on

 plot([0,max(residue)],SHIFT_med,'r')
 plot([0,max(residue)],SHIFT_med+SHIFT_std,'r-.')
 plot([0,max(residue)],SHIFT_med-SHIFT_std,'r-.')
 plot([0,max(residue)],SHIFT_med+2*SHIFT_std,'r:')

text(max(residue)/3,max(SHIFT_purge)*0.8,num2str(SHIFT_med(1)))
text(max(residue)/4,max(SHIFT_purge)*0.8,'med')
text(max(residue)/3,max(SHIFT_purge)*0.6,num2str(SHIFT_std))
text(max(residue)/4,max(SHIFT_purge)*0.6,'std')
text(max(residue)/3,max(SHIFT_purge)*0.9,num2str(max(SHIFT_purge)))
text(max(residue)/4,max(SHIFT_purge)*0.9,'max')
figure(m+100)
clf
bar(residue_purge,SHIFTH_purge)
title('only H')
figure(m+200)
clf
bar(residue_purge,SHIFTN_purge)
title('only N')
%}

 
for m = 2:n      %output loop after purge
    
figure(300)     %sub plot for summary HN
SHIFT_med(m)=median(SHIFT_purge(:,m));
SHIFT_mean(m)=mean(SHIFT_purge(:,m));
SHIFT_std(m)=std(SHIFT_purge(:,m));
subplot(nc,3,m)
bar(residue_purge,SHIFT_purge(:,m))
axis tight
hold on
plot([0,max(residue)],[SHIFT_med(m) SHIFT_med(m)],'r')
plot([0,max(residue)],[(SHIFT_med(m)+SHIFT_std(m)) (SHIFT_med(m)+SHIFT_std(m))],'r-.')
plot([0,max(residue)],[(SHIFT_med(m)-SHIFT_std(m)) (SHIFT_med(m)-SHIFT_std(m))],'r-.')  
plot([0,max(residue)],[(SHIFT_med(m)+2*SHIFT_std(m)) (SHIFT_med(m)+2*SHIFT_std(m))],'r:')
titre=['HN' num2str(percentage(m), '%10.1f\n') '%'];
title(titre)
text(max(residue)*0.5,max(SHIFT_purge(:,m))*0.8,num2str(SHIFT_med(m)))
text(max(residue)*0.1,max(SHIFT_purge(:,m))*0.8,'med')
text(max(residue)*0.5,max(SHIFT_purge(:,m))*0.6,num2str(SHIFT_std(m)))
text(max(residue)*0.1,max(SHIFT_purge(:,m))*0.6,'std')
text(max(residue)*0.5,max(SHIFT_purge(:,m))*0.9,num2str(max(SHIFT_purge(:,m))))
text(max(residue)*0.1,max(SHIFT_purge(:,m))*0.9,'max')

%plot for current titration point
 figure(m)
 clf
bar(residue_purge,SHIFT_purge(:,m))
axis tight
Out1=find(SHIFT_purge(:,m)>(SHIFT_med(m)+SHIFT_std(m)));
Out2=find(SHIFT_purge(:,m)>(SHIFT_med(m)+2*SHIFT_std(m)));
hold on
bar(residue_purge(Out1),SHIFT_purge(Out1,m),'r')
% bar(residue_purge(Out2),SHIFT_purge(Out2,m),'y')


xlim(resrange);
if (shiftY)
    ylim([0 shiftY])
end

hold on
plot([0,max(residue)],[SHIFT_med(m) SHIFT_med(m)],'r')
plot([0,max(residue)],[(SHIFT_med(m)+SHIFT_std(m)) (SHIFT_med(m)+SHIFT_std(m))],'r-.')
plot([0,max(residue)],[(SHIFT_med(m)-SHIFT_std(m)) (SHIFT_med(m)-SHIFT_std(m))],'r-.')  
plot([0,max(residue)],[(SHIFT_med(m)+2*SHIFT_std(m)) (SHIFT_med(m)+2*SHIFT_std(m))],'r:')
titre=['HN' num2str(percentage(m), '%10.1f\n') '%'];
title(titre)
text(max(residue)*0.5,max(SHIFT_purge(:,m))*0.8,num2str(SHIFT_med(m)))
text(max(residue)*0.1,max(SHIFT_purge(:,m))*0.8,'med')
text(max(residue)*0.5,max(SHIFT_purge(:,m))*0.6,num2str(SHIFT_std(m)))
text(max(residue)*0.1,max(SHIFT_purge(:,m))*0.6,'std')
text(max(residue)*0.5,max(SHIFT_purge(:,m))*0.9,num2str(max(SHIFT_purge(:,m))))
text(max(residue)*0.1,max(SHIFT_purge(:,m))*0.9,'max')

% figure(301)     %sub plot for summary N
% subplot(nc,3,m)
% bar(residue_purge,SHIFTN_purge)
% axis tight
% titre=['N' num2str(PC, '%10.1f\n') '%'];
% title(titre)
% figure(302)     %sub plot for summary H
% subplot(nc,3,m)
% bar(residue_purge,SHIFTH_purge)
% axis tight
% titre=['H' num2str(PC, '%10.1f\n') '%'];
% title(titre)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%BREAKKKKKKKKKKKKKKKKKKKKK coding
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%06/20/16%%%%%%%%%%%%%%%%%%%


%save contacts only based on H/N
contact = residue_purge(find(SHIFT_purge(:,m)>(SHIFT_med(m)+SHIFT_std(m))));
SHIFTcontact=SHIFT_purge(find(SHIFT_purge(:,m)>(SHIFT_med(m)+SHIFT_std(m))));
contact2 = residue_purge(find(SHIFT_purge(:,m)>(SHIFT_med(m)+2*SHIFT_std(m))));
SHIFTcontact2 = SHIFT_purge(find(SHIFT_purge(:,m)>(SHIFT_med(m)+2*SHIFT_std(m))));
savefile=['contact' num2str(percentage(m), '%10.1f\n')];
savefile2=['contact2_' num2str(percentage(m), '%10.1f\n')];
save(savefile, 'contact', '-ASCII')
save(savefile2, 'contact2', '-ASCII')
%print to file figures for H, N and H/N
figurename=['overall' num2str(percentage(m), '%10.1f\n') '.pdf'];
fhandel=['-f' int2str(m)];
print(fhandel, '-dpdf', figurename)
%figurename=['overallH' num2str(PC, '%10.1f\n') '.pdf'];
%fhandel=['-f' int2str(m+100)];
%print(fhandel, '-dpdf', figurename)
%figurename=['overallN' num2str(PC, '%10.1f\n') '.pdf'];
%fhandel=['-f' num2str(m+200)];
%print(fhandel, '-dpdf', figurename)
%store all shifts to be ploted on a residue per residue basis
% ALLSHIFT=[ALLSHIFT SHIFT_purge];
% percentage=[percentage PC];

%export all residues except bad residues (overlap etc.) use as input for
%replaceBfactor
res_shift=[residue_purge SHIFT_purge(:,m)];
filename=['shift_HN_purge_' num2str(percentage(m), '%10.1f\n')];
fid=fopen(filename, 'w');
fprintf(fid,'%3u %11.9f\n',res_shift');
fclose(fid);
%export contact residues for 1std use as input for
%replaceBfactor
res_shift=[contact SHIFTcontact];
filename=['shift_HN_contact_' num2str(percentage(m), '%10.1f\n')];
fid=fopen(filename, 'w');
fprintf(fid,'%3u %11.9f\n',res_shift');
fclose(fid);
%export contact residues for 1std use as input for
%replaceBfactor
res_shift=[contact2 SHIFTcontact2];
filename=['shift_HN_contact2_' num2str(percentage(m), '%10.1f\n')];
fid=fopen(filename, 'w');
fprintf(fid,'%3u %11.9f\n',res_shift');
fclose(fid);

end %output loop

figure(300)
print -dpdf overall_summaryHN.pdf
%figure(301)
%print -dpdf overall_summaryN.pdf
%figure(302)
%print -dpdf overall_summaryH.pdf
%close all
i=0;
b=1;
fig_count=0;
figure(1+1000)
%begin analysis on a per residue basis, shifts
for i=1:size(INT,1)
    %figure(i+1000)
    subplot(5,3,b)
    plot(percentage,SHIFT_purge(i,:),'gd-')
     hold on
    %fit intial point linear regime
    [fit_shift, g_shift]=fit(percentage(1:3)',SHIFT_purge(i,1:3)', 'poly1');
    shift_fit(i,1)=fit_shift.p1;
    shift_fit(i,2)=fit_shift.p2;
    shift_fit(i,3)=g_shift.rsquare;
    
    %plot(percentage(1:3), (percentage(1:3).*fit_shift.p1+fit_shift.p2),'b-');
    titre=['r' Sample int2str(residue_purge(i))];
    title(titre)
    %axis auto
    %axis([min(percentage) max(percentage) 0 max(max(SHIFT_purge))])
    axis([min(percentage) max(percentage) 0 0.14])
    
%     subplot(2,3,b+3)
%     plot(percentage,SHIFT_purge(i,:),'gd-')
%     titre=['sc. r ' int2str(residue_purge(i))];
%     title(titre)
%     axis([min(percentage) max(percentage) 0 max(max(SHIFT_purge))])
%     i;
%     b;
    if b == 15
        figurename=['SHIFTresidue_fit' int2str(residue_purge(i-b+1)) 'to' int2str(residue_purge(i))];
        print('-dpdf',figurename,'-fillpage')
        %clf
        close
        b=0;
        fig_count=fig_count+1;  %for debugging
        figure(i+1000+1)
    end
    b=b+1;
end
%  a = axes;
%  t1 = title(Sample);
%  a.Visible = 'off'; % set(a,'Visible','off');
%  t1.Visible = 'on'; % set(t1,'Visible','on');

figurename=['SHIFTresidue_fit' int2str(residue(i-b+1)) 'to' int2str(residue(i))];
print('-dpdf',figurename)
%clf


% i=0;
% b=1;
% fig_count=0;
% figure(1+2000)
% %begin analysis on a per residue basis, INT
% for i=1:size(INT,1)
%     %figure(i+1000)
%     subplot(2,3,b)
%     %plot(percentage,D_INT(i,:),'gd-')
%     plot(percentage,DN_INT(i,:),'gd-')
%      hold on
%     %fit intial point linear regime
%     [fit_int, g_int]=fit(percentage(1:3)',DN_INT(i,1:3)', 'poly1');
%     int_fit(i,1)=fit_int.p1;
%     int_fit(i,2)=fit_int.p2;
%     int_fit(i,3)=g_int.rsquare;
%     plot(percentage(1:3), (percentage(1:3).*fit_int.p1+fit_int.p2),'b-');
%     titre=['DN r ' int2str(residue_purge(i)) ' chi2 =' int2str(g_int.rsquare)];
%     title(titre)
%     axis([min(percentage) max(percentage) min(min(D_INT)) max(max(D_INT))])
    %axis auto
    
%     subplot(2,3,b+3)
%     plot(percentage,N_INT(i,:),'gd-')
%     titre=['I/Io r ' int2str(residue_purge(i))];
%     title(titre)
%     axis([min(percentage) max(percentage) min(min(N_INT)) max(max(N_INT))])
%     %axis auto
%     i;
%     b;
%     if b == 3
%         figurename=['DN_INTresidue_fit' int2str(residue_purge(i-b+1)) 'to' int2str(residue_purge(i))];
%         %figurename=['INTresidue' int2str(residue_purge(i-b+1)) 'to' int2str(residue_purge(i))];
%         print('-dpdf',figurename)
%         %clf
%         close
%         b=0;
%         fig_count=fig_count+1;  %for debugging
%         figure(i+2000+1)
%     end
%     b=b+1;
% end
% figurename=['DN_INTresidue_fit' int2str(residue(i-b+1)) 'to' int2str(residue(i))];
% print('-dpdf',figurename)


% a1=[18 32]-res_off;
% l1=[33 53]-res_off;
% a2=[54 66]-res_off;
% l2=[67 72]-res_off;
% a3=[73 78]-res_off;
% l3=[79 81]-res_off;
% a4=[82 89]-res_off;

%Here you must account for offset; can be fixed by using vectors above and
%modify definitions below.

if(0)

a1=residue(find(residue>17&residue<33));
a2=residue(find(residue>53&residue<67));
a3=residue(find(residue>72&residue<79));
a4=residue(find(residue>81&residue<90));
l1=residue(find(residue>32&residue<54));
l2=residue(find(residue>66&residue<73));
l3=residue(find(residue>78&residue<82));

shift_a1=shift_fit((find(residue>17&residue<33)),1);
shift_a2=shift_fit((find(residue>53&residue<67)),1);
shift_a3=shift_fit((find(residue>72&residue<79)),1);
shift_a4=shift_fit((find(residue>81&residue<90)),1);
shift_l1=shift_fit((find(residue>32&residue<54)),1);
shift_l2=shift_fit((find(residue>66&residue<73)),1);
shift_l3=shift_fit((find(residue>78&residue<82)),1);

int_a1=int_fit((find(residue>17&residue<33)),1);
int_a2=int_fit((find(residue>53&residue<67)),1);
int_a3=int_fit((find(residue>72&residue<79)),1);
int_a4=int_fit((find(residue>81&residue<90)),1);
int_l1=int_fit((find(residue>32&residue<54)),1);
int_l2=int_fit((find(residue>66&residue<73)),1);
int_l3=int_fit((find(residue>78&residue<82)),1);




figure(4001)     %sub plot for summary HN
subplot(2,1,1)
shift_fit_med=median(shift_fit(:,1));
shift_fit_mean=mean(shift_fit(:,1));
shift_fit_std=std(shift_fit(:,1));
%bar(residue_purge,shift_fit(:,1))
xlim(resrange);
if (shiftslopeY)
    ylim([0 shiftslopeY])
end

bar(residue_purge,shift_fit(:,1),'w')
hold on
bar(a1,shift_a1,0.8,'b')
hold on
bar(a2,shift_a2,'g')
hold on
bar(a3,shift_a3,'m')
hold on
bar(a4,shift_a4,'c')
hold on
bar(l1,shift_l1,'r')
hold on
bar(l2,shift_l2,'k')
hold on
%bar(l3(1),shift_l3(1),0.8,'y')  %matlab bug for width of 2 entries (1, or >3 are fine)
bar(l3(2),shift_l3(2),0.8,'y')
hold on

%axis tight

hold on
plot([0,max(residue)],[shift_fit_med shift_fit_med],'r')
plot([0,max(residue)],[(shift_fit_med+shift_fit_std) (shift_fit_med+shift_fit_std)],'r-.')
plot([0,max(residue)],[(shift_fit_med-shift_fit_std) (shift_fit_med-shift_fit_std)],'r-.')  
plot([0,max(residue)],[(shift_fit_med+2*shift_fit_std) (shift_fit_med+2*shift_fit_std)],'r:')
%titre=['CSP slope' num2str(percentage(m), '%10.1f\n') '%'];
title('CSP slope')
text(max(residue)*0.5,max(shift_fit(:,1))*0.8,num2str(shift_fit_med))
text(max(residue)*0.1,max(shift_fit(:,1))*0.8,'med')
text(max(residue)*0.5,max(shift_fit(:,1))*0.6,num2str(shift_fit_std))
text(max(residue)*0.1,max(shift_fit(:,1))*0.6,'std')
text(max(residue)*0.5,max(shift_fit(:,1))*0.9,num2str(max(shift_fit(:,1))))
text(max(residue)*0.1,max(shift_fit(:,1))*0.9,'max')
%axis tight
%xlim(resrange);
%axes;

subplot(2,1,2)
int_fit_med=median(int_fit(:,1));
int_fit_mean=mean(int_fit(:,1));
int_fit_std=std(int_fit(:,1));
%bar(residue_purge,int_fit(:,1)) %comment this and insert SS-based color code below.
xlim(resrange);
if (intslopeY)
    ylim([0 intslopeY])
end
bar(residue_purge,int_fit(:,1),'w')

hold on
bar(a1,int_a1,'b')
hold on
bar(a2,int_a2,'g')
hold on
bar(a3,int_a3,'m')
hold on
bar(a4,int_a4,'c')
hold on
bar(l1,int_l1,'r')
hold on
bar(l2,int_l2,'k')
hold on
bar(l3(1),int_l3(1),'y')
bar(l3(2),int_l3(2),'y')
hold on
%axis tight
%xlim(resrange);
hold on
plot([0,max(residue)],[int_fit_med int_fit_med],'r')
plot([0,max(residue)],[(int_fit_med+int_fit_std) (int_fit_med+int_fit_std)],'r-.')
plot([0,max(residue)],[(int_fit_med-int_fit_std) (int_fit_med-int_fit_std)],'r-.')  
plot([0,max(residue)],[(int_fit_med+2*int_fit_std) (int_fit_med+2*int_fit_std)],'r:')
%titre=['HN' num2str(percentage(m), '%10.1f\n') '%'];
title('(Io-Ix)/Io slope')
text(max(residue)*0.5,max(int_fit(:,1))*0.8,num2str(int_fit_med))
text(max(residue)*0.1,max(int_fit(:,1))*0.8,'med')
text(max(residue)*0.5,max(int_fit(:,1))*0.6,num2str(int_fit_std))
text(max(residue)*0.1,max(int_fit(:,1))*0.6,'std')
text(max(residue)*0.5,max(int_fit(:,1))*0.9,num2str(max(int_fit(:,1))))
text(max(residue)*0.1,max(int_fit(:,1))*0.9,'max')
%axis tight
%xlim(resrange);
figurename='CSP_INT_stacked';
print('-dpdf',figurename);


figure(4000)
plot(shift_fit(:,1),int_fit(:,1),'wo','MarkerFaceColor','w','MarkerEdgeColor','k')
hold on
%plot(shift_a1(:,1),int_a1(:,1),'bo')
plot(shift_a1(:,1),int_a1(:,1),'bo','MarkerFaceColor','b')
hold on
plot(shift_a2(:,1),int_a2(:,1),'go','MarkerFaceColor','g')
hold on
plot(shift_a3(:,1),int_a3(:,1),'mo','MarkerFaceColor','m')
hold on
plot(shift_a4(:,1),int_a4(:,1),'co','MarkerFaceColor','c')
hold on

plot(shift_l1(:,1),int_l1(:,1),'or','MarkerFaceColor','r')
hold on
plot(shift_l2(:,1),int_l2(:,1),'ok','MarkerFaceColor','k')
hold on
plot(shift_l3(:,1),int_l3(:,1),'oy','MarkerFaceColor','y')
hold on


plot([0,max(shift_fit(:,1))],[int_fit_med int_fit_med],'r')
plot([0,max(shift_fit(:,1))],[(int_fit_med+int_fit_std) (int_fit_med+int_fit_std)],'r-.')
plot([shift_fit_med shift_fit_med],[0,max(int_fit(:,1))],'r')
plot([(shift_fit_med+shift_fit_std) (shift_fit_med+shift_fit_std)],[0,max(int_fit(:,1))],'r-.')

resstr=num2str(residue_purge);
cellres=cellstr(resstr);
text(shift_fit(:,1),int_fit(:,1),cellres)
title('shift/int variation');
axis auto
figurename='CSP_INT_corr';
print('-dpdf',figurename);

end

%clf

%syms Ka Ptot Ltot L
%P=Ka*L^2+(1+Ka*(Ptot-Ltot))*L-Ltot
%Li=solve(P,L)
%syms delta detla_p detla_pl
%delta=(detla_p + detla_pl*Ka*Li(1))/(1+Ka*Li(1))

%close
 
