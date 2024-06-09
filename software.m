clc;
clear all;
close all

%%%%%Es季节变化
load('xamm.mat');load('xxmm.mat');load('wlmqmm.mat');load('szmm.mat');load('qdmm.mat');
load('mzlmm.mat');load('lzmm.mat');load('lsmm.mat');load('kmmm.mat');load('hkmm.mat');
load('gzmm.mat');load('bjmm.mat');load('ccmm.mat');load('cqmm.mat');load('okimm.mat');
load('kokmm.mat');load('akimm.mat');load('yammm.mat');load('wakmm.mat');load('ssmm.mat');
load('whmm.mat');
for i=1:12
    wl(i,:)=nanmean(wlmqmm(i:12:length(wlmqmm),:));
    bj(i,:)=nanmean(bjmm(i:12:length(bjmm),:));
    mz(i,:)=nanmean(mzlmm(i:12:length(mzlmm),:));
    cc(i,:)=nanmean(ccmm(i:12:length(ccmm),:));
    wak(i,:)=nanmean(wakmm(i:12:length(wakmm),:));
    lz(i,:)=nanmean(lzmm(i:12:length(lzmm),:));
    xx(i,:)=nanmean(xxmm(i:12:length(xxmm),:));
    qd(i,:)=nanmean(qdmm(i:12:length(qdmm),:));
    aki(i,:)=nanmean(akimm(i:12:length(akimm),:));
    kok(i,:)=nanmean(kokmm(i:12:length(kokmm),:));
    ls(i,:)=nanmean(lsmm(i:12:length(lsmm),:));
    cq(i,:)=nanmean(cqmm(i:12:length(cqmm),:));
    wh(i,:)=nanmean(whmm(i:12:length(whmm),:));
    ss(i,:)=nanmean(ssmm(i:12:length(ssmm),:));
    yam(i,:)=nanmean(yammm(i:12:length(yammm),:));
    km(i,:)=nanmean(kmmm(i:12:length(kmmm),:));
    gz(i,:)=nanmean(gzmm(i:12:length(gzmm),:));
    hk(i,:)=nanmean(hkmm(i:12:length(hkmm),:));
    sz(i,:)=nanmean(szmm(i:12:length(szmm),:));
    oki(i,:)=nanmean(okimm(i:12:length(okimm),:));
    xa(i,:)=nanmean(xamm(i:12:length(xamm),:));
end
bjs=nanmean(nanmax(bj(5:7,:)));
ccs=nanmean(nanmax(cc(5:7,:)));
cqs=nanmean(nanmax(cq(5:7,:)));
hks=nanmean(nanmax(hk(5:7,:)));
wls=nanmean(nanmax(wl(5:7,:)));
lss=nanmean(nanmax(ls(5:7,:)));
lzs=nanmean(nanmax(lz(5:7,:)));
sss=nanmean(nanmax(ss(5:7,:)));
xxs=nanmean(nanmax(xx(5:7,:)));
mzs=nanmean(nanmax(mz(5:7,:)));
szs=nanmean(nanmax(sz(5:7,:)));
whs=nanmean(nanmax(wh(5:7,:)));
qds=nanmean(nanmax(qd(5:7,:)));
gzs=nanmean(nanmax(gz(5:7,:)));
kms=nanmean(nanmax(km(5:7,:)));
waks=nanmean(nanmax(wak(5:7,:)));
yams=nanmean(nanmax(yam(5:7,:)));
akis=nanmean(nanmax(aki(5:7,:)));
koks=nanmean(nanmax(kok(5:7,:)));
okis=nanmean(nanmax(oki(5:7,:)));
xas=nanmean(nanmax(xa(5:7,:)));

bjw(1,:)=bj(1,:);
bjw(2:3,:)=bj(11:12,:);
bjw=nanmean(nanmean(bjw));

ccw(1,:)=cc(1,:);
ccw(2:3,:)=cc(11:12,:);
ccw=nanmean(nanmin(ccw));

cqw(1,:)=cq(1,:);
cqw(2:3,:)=cq(11:12,:);
cqw=nanmean(nanmin(cqw));

hkw(1,:)=hk(1,:);
hkw(2:3,:)=hk(11:12,:);
hkw=nanmean(nanmin(hkw));

wlw(1,:)=wl(1,:);
wlw(2:3,:)=wl(11:12,:);
wlw=nanmean(nanmin(wlw));

lsw(1,:)=ls(1,:);
lsw(2:3,:)=ls(11:12,:);
lsw=nanmean(nanmin(lsw));

lzw(1,:)=lz(1,:);
lzw(2:3,:)=lz(11:12,:);
lzw=nanmean(nanmin(lzw));

ssw(1,:)=ss(1,:);
ssw(2:3,:)=ss(11:12,:);
ssw=nanmean(nanmin(ssw));

mzw(1,:)=mz(1,:);
mzw(2:3,:)=mz(11:12,:);
mzw=nanmean(nanmin(mzw));

szw(1,:)=sz(1,:);
szw(2:3,:)=sz(11:12,:);
szw=nanmean(nanmin(szw));

whw(1,:)=wh(1,:);
whw(2:3,:)=wh(11:12,:);
whw=nanmean(nanmin(whw));

qdw(1,:)=qd(1,:);
qdw(2:3,:)=qd(11:12,:);
qdw=nanmean(nanmin(qdw));

gzw(1,:)=gz(1,:);
gzw(2:3,:)=gz(11:12,:);
gzw=nanmean(nanmin(gzw));

xxw(1,:)=xx(1,:);
xxw(2:3,:)=xx(11:12,:);
xxw=nanmean(nanmin(xxw));

kmw(1,:)=km(1,:);
kmw(2:3,:)=km(11:12,:);
kmw=nanmean(nanmin(kmw));

wakw(1,:)=wak(1,:);
wakw(2:3,:)=wak(11:12,:);
wakw=nanmean(nanmin(wakw));

yamw(1,:)=yam(1,:);
yamw(2:3,:)=yam(11:12,:);
yamw=nanmean(nanmin(yamw));

akiw(1,:)=aki(1,:);
akiw(2:3,:)=aki(11:12,:);
akiw=nanmean(nanmin(akiw));

kokw(1,:)=kok(1,:);
kokw(2:3,:)=kok(11:12,:);
kokw=nanmean(nanmin(kokw));

okiw(1,:)=oki(1,:);
okiw(2:3,:)=oki(11:12,:);
okiw=nanmean(nanmin(okiw));

xaw(1,:)=xa(1,:);
xaw(2:3,:)=xa(11:12,:);
xaw=nanmean(nanmin(xaw));
%%%%%Es季节变化SSS

load('bjtl.mat');load('bjtu.mat');load('whtl.mat');load('whtu.mat');
load('wlmqtl.mat');load('wlmqtu.mat');load('lstl.mat');load('lstu.mat');
load('lztl.mat');load('lztu.mat');load('xatl.mat');load('xatu.mat');
load('xxtl.mat');load('xxtu.mat');load('qdtl.mat');load('qdtu.mat');
load('mzltl.mat');load('mzltu.mat');load('cctl.mat');load('cctu.mat');
load('cqtl.mat');load('cqtu.mat');load('sztl.mat');load('sztu.mat');
load('kmtl.mat');load('kmtu.mat');load('gztl.mat');load('gztu.mat');
load('hktl.mat');load('hktu.mat');load('akitl.mat');load('akitu.mat');
load('okitl.mat');load('okitu.mat');load('yamtl.mat');load('yamtu.mat');
load('waktl.mat');load('waktu.mat');load('koktl.mat');load('koktu.mat');
bjt=(bjtl+bjtu)/2;wht=(whtl+whtu)/2;wlmqt=(wlmqtl+wlmqtu)/2;lst=(lstl+lstu)/2;
lzt=(lztl+lztu)/2;xat=(xatl+xatu)/2;xxt=(xxtl+xxtu)/2;qdt=(qdtl+qdtu)/2;
mzlt=(mzltl+mzltu)/2;cct=(cctl+cctu)/2;cqt=(cqtl+cqtu)/2;szt=(sztl+sztu)/2;
kmt=(kmtl+kmtu)/2;gzt=(gztl+gztu)/2;hkt=(hktl+hktu)/2;akit=(akitl+akitu)/2;
okit=(okitl+okitu)/2;yamt=(yamtl+yamtu)/2;wakt=(waktl+waktu)/2;kokt=(koktl+koktu)/2;
bjtm=max(bjt);whtm=max(wht);wlmqtm=max(wlmqt);lstm=max(lst);
lztm=max(lzt);xatm=max(xat);xxtm=max(xxt);qdtm=max(qdt);
mzltm=max(mzlt);cctm=max(cct);cqtm=max(cqt);sztm=max(szt);
kmtm=max(kmt);gztm=max(gzt);hktm=max(hkt);akitm=max(akit);
okitm=max(okit);yamtm=max(yamt);waktm=max(wakt);koktm=max(kokt);

lonB=70:.5:150;
latB=15:.5:55;

latA=[49.58  43.84  43.75  40.11  36.24  36.06  31     29.63 29.5  25.5  23.15   20    30.5  31.2  26.3  39.7  45.4 35.26  31.3  35.7 34.23];
lonA=[117.45 125.27 87.63 116.27 120.41 103.87 121.24 91.28 106.4 103.8 113.35 110.33 114.4 130.6 127.8 140.1 141.7 113.9  120.6 139.5 108.92];
Zday=[mzltm cctm wlmqtm bjtm qdtm lztm sztm lstm cqtm kmtm gztm hktm whtm yamtm okitm akitm waktm  xxtm sztm koktm xatm]; 
Zw=[mzs ccs wls bjs qds lzs sss lss cqs kms gzs hks whs yams okis akis waks  xxs szs koks xas]; 
Z=zeros(81,161);
for i=1:81
    for j=1:161
        count=0;
        for ii=1:21
            d(ii)=((lonA(ii)-lonB(j))^2+25*(latA(ii)-latB(i))^2)^1/2;
            Z(i,j)=Z(i,j)+1/d(ii)*Zday(ii);
            count=count+1/d(ii);
        end
        Z(i,j)=Z(i,j)/count;
    end
end
Z1=Z;
Z=zeros(81,161);
for i=1:81
    for j=1:161
        count=0;
        for ii=1:21
            d(ii)=((lonA(ii)-lonB(j))^2+25*(latA(ii)-latB(i))^2)^1/2;
            Z(i,j)=Z(i,j)+1/d(ii)*Zw(ii);
            count=count+1/d(ii);
        end
        Z(i,j)=Z(i,j)/count;
    end
end
Z2=Z;
% v=0.2:.05:.7;
% contourf(70:140,15:55,Z,v);
%  Z1=Z;
figure
subplot(2,1,2), 
% contourf(70:.5:150,15:.5:55,Z1,96,'linestyle','none')
imagesc(70:.5:150,15:.5:55,Z1)
colorbar;
set(gca,'ydir','normal')
grid;
ylabel('Latitude(deg.)'); xlabel('Longitude(deg.)');hold on
% plot_map;
%%%%%%画地图
Qingdao=[36.24,120.41];
Suzhou=[31.32,120.63];
Xinxiang=[35.26,113.95];
Beijing=[40.11,116.27];
Changchun=[43.84,125.27];
Manzhouli=[49.58,117.45];
Lanzhou=[36.06,103.87];
Chongqing=[29.5,106.4];
Kunming=[25.5,103.8];
Guangzhou=[23.15,113.35];
Haikou=[20,110.33];
Lhasa=[29.63,91.28];
Urumchi=[43.75,87.63];
Wuhan=[30.5,114.4];
Xian=[34.23,108.92];
Akita=[39.7,140.10];
Okinawa=[26.3,127.80];
Yamagawa=[31.2,130.60];
Wakkanai=[45.4,141.7];
Koku=[35.7,139.5];

load('worldMaps.mat')
longitude=ll_world(:,1);
latitude=ll_world(:,2);
plot(ll_world(:,2),ll_world(:,1),'k','LineWidth',0.5);
if(1)
% figure,
hold on,plot(longitude,latitude,'b.','MarkerSize',.001)
hold on,plot(Qingdao(1,2),Qingdao(1,1),'ro','MarkerSize',6,'MarkerFaceColor','r'),
hold on,plot(Suzhou(1,2),Suzhou(1,1),'ro','MarkerSize',6,'MarkerFaceColor','r'),
hold on,plot(Beijing(1,2),Beijing(1,1),'ro','MarkerSize',6,'MarkerFaceColor','r'),
hold on,plot(Xinxiang(1,2),Xinxiang(1,1),'ro','MarkerSize',6,'MarkerFaceColor','r'),
hold on,plot(Changchun(1,2),Changchun(1,1),'ro','MarkerSize',6,'MarkerFaceColor','r'),
hold on,plot(Manzhouli(1,2),Manzhouli(1,1),'ro','MarkerSize',6,'MarkerFaceColor','r'),
hold on,plot(Lanzhou(1,2),Lanzhou(1,1),'ro','MarkerSize',6,'MarkerFaceColor','r'),
hold on,plot(Chongqing(1,2),Chongqing(1,1),'ro','MarkerSize',6,'MarkerFaceColor','r'),
hold on,plot(Kunming(1,2),Kunming(1,1),'ro','MarkerSize',6,'MarkerFaceColor','r'),
hold on,plot(Haikou(1,2),Haikou(1,1),'ro','MarkerSize',6,'MarkerFaceColor','r'),
hold on,plot(Lhasa(1,2),Lhasa(1,1),'ro','MarkerSize',6,'MarkerFaceColor','r'),
hold on,plot(Urumchi(1,2),Urumchi(1,1),'ro','MarkerSize',6,'MarkerFaceColor','r'),
hold on,plot(Guangzhou(1,2),Guangzhou(1,1),'ro','MarkerSize',6,'MarkerFaceColor','r'),
hold on,plot(Wuhan(1,2),Wuhan(1,1),'ro','MarkerSize',6,'MarkerFaceColor','r'),
hold on,plot(Xian(1,2),Xian(1,1),'ro','MarkerSize',6,'MarkerFaceColor','r'),
hold on,plot(Akita(1,2),Akita(1,1),'ro','MarkerSize',6,'MarkerFaceColor','r'),
hold on,plot(Okinawa(1,2),Okinawa(1,1),'ro','MarkerSize',6,'MarkerFaceColor','r'),
hold on,plot(Yamagawa(1,2),Yamagawa(1,1),'ro','MarkerSize',6,'MarkerFaceColor','r'),
hold on,plot(Wakkanai(1,2),Wakkanai(1,1),'ro','MarkerSize',6,'MarkerFaceColor','r'),
hold on,plot(Koku(1,2),Koku(1,1),'ro','MarkerSize',6,'MarkerFaceColor','r'),

set(gca,'ylim',[15,55],'xlim',[70,150]);
set(gca,'xtick',70:10:150,'xticklabel',{'70','80','90','100','110','120','130','140','150'}); 
set(gca,'ytick',15:10:55,'yticklabel',{'15','25','35','45','55'}); 
xlabel('Longitude/deg.'); ylabel('Latitude/deg.');

grid
end
set(gca,'ydir','normal')

%%%%%%画地图SSS
subplot(2,1,1), 
imagesc(70:.5:150,15:.5:55,Z2)
colorbar;
grid;
set(gca,'ydir','normal')
ylabel('Latitude(deg.)'); xlabel('Longitude(deg.)');hold on
% plot_map;
% save('Z1.mat','Z1');

% clear bjw; clear bjs
%  %%%%%%%%%%%%%%%%%%%%%季节变化，暂时变灰，求季节变化时重新启用（终止符）


