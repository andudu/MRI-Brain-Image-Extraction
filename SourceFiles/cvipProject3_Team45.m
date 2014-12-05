function cvipProject3_Team45()

%CVIP Project 3
%By Priyanka Kulkarni
%Gaurav Talokar
%Team 45
    
%----------Edge Preserve Smoothing-------------------------%

I=imread('test.bmp');
Img=im2double(I);

figure(1),imshow(I),title('original');

[gx, gy] = imgradientxy(Img);
[Gmag, Gdir] = imgradient(gx, gy);      %Compute the gradient magnitude
k = 5;
win=exp(-1*abs(Gmag.^2)/(2*k*k));
               
[row,col]=size(Img);        %Get Image Size

for iterations=1:10
 for xx=2:row-1         %for number of rows     
     for yy=2:col-1     %For number of columns
         
            N_val=0;
            TempI=0;
            
             %Summation Computation
              for i=-1 :1
                for j=-1 :1
                    N_val=N_val+ win(xx+i,yy+j);  %computations 
                    TempI=TempI + Img(xx+i, yy+j)*win(xx+i,yy+j);
                 end
              end
                EdgeP_Image(xx,yy)=TempI/N_val;
      end
 end
end
 Img=EdgeP_Image;
 figure(2);
imshow(Img), title('Edge Preserve Smoothing');


%----------Canny Edge Detection---------------------------
t=0.04;
Canny=edge(I,'canny',t);
figure(3);
imshow(Canny), title('Canny Edge Image');

%--------Region Growing Using 4-neighbourhood--------------
%Ig=imread('test.bmp');
[row,col]=size(Img);
NewI=Img;

thresh=0.0035;
Z=zeros(256);



val=1;  
R=zeros(256);

count=1;
c=1;
for i=2:row-1
    
    for j=2:col-1
             
             cur_val=NewI(i,j);
             neighbour_val2=NewI(i,j-1);       %left
             neighbour_val3=NewI(i,j+1);       %Right
             neighbour_val1=NewI(i+1,j);       %down
             neighbour_val4=NewI(i-1,j);       %up

             if((cur_val-neighbour_val2<thresh)&&(cur_val-neighbour_val3 <thresh)&&(cur_val-neighbour_val4<thresh)&&(cur_val-neighbour_val1<thresh))
                 Z(i,j-1)=0;
                 
                     Z(i,j+1)= NewI(i,j);
                     Z(i-1,j)=NewI(i,j);
                     Z(i+1,j)=NewI(i,j);
                     Z(i,j)=NewI(i,j);
                     Z(i,j-1)=NewI(i,j);
                     
                     
                  %checking if LEFT pixel is lbeled or not
                    if(R(i,j-1)>0)
                       R(i,j)=R(i,j-1);     %center
                        R(i,j+1)=R(i,j-1);  %right
                       R(i-1,j)=R(i,j-1);  %up
                       R(i+1,j)=R(i,j-1);  %down
                      
                    end
                   
                %checking if RIGHT pixel is lbeled or not
                   if(R(i,j+1)>0)
                       R(i,j)=R(i,j+1);   %center
                       R(i,j-1)=R(i,j+1);  %left
                       R(i-1,j)=R(i,j+1);   %up
                       R(i+1,j)=R(i,j+1);   %down
                      
                   end
                   
               %checking if DOWNWARD pixel is lbeled or not
                   if(R(i+1,j)>0)
                       
                       R(i,j)=R(i+1,j);    %center
                       R(i,j-1)=R(i+1,j);   %left
                       R(i,j+1)=R(i+1,j);   %Right
                       R(i-1,j)=R(i+1,j);    %up
                      
                   end
                   
              %checking if UPPER pixel is lbeled or not
                  if(R(i-1,j)>0)    
                       R(i,j)=R(i-1,j);     %center
                       R(i,j-1)=R(i-1,j);    %left
                       R(i,j+1)=R(i-1,j);    %Right
                       R(i+1,j)=R(i-1,j);     %down
                       
                   else
                      val=val+1;
                      R(i,j)=val;
                      R(i-1,j)=val; 
                      R(i+1,j)=val; 
                      R(i,j-1)=val; 
                      R(i,j+1)=val; 
                     
               end
                                   
         end           
    end
    

end


figure(4),imshow(Z),title('region growing');


%--------------Region Merging---------------------------------%
total_reg=val;
arr=zeros(1,256*256);
%To store standard Deviations
regionstdeviation=zeros(1,total_reg);
%To compute the mean of all region intensisties
regionmeans=zeros(1,total_reg);
%To number of regions sharing labels
regioncount=zeros(1,total_reg);
%TO store region length
length=zeros(1,total_reg);


 for i=1:total_reg
    count=1;
    cost=length(1,i);
    for ii=2:row-1
         for jj=2:col-1
               reg=R(ii,jj); 
              if(reg==i)
                 arr(1,count)=Z(ii,jj);
                 count=count+1;
              end

         end
    end
    
    
    arraytrim=arr(1:count);
    regionmeans(1,i)=mean(arraytrim);
    regioncount(1,i)=count;
    regionstdeviation(1,i)=std(arraytrim);
 end
 
 countregion=zeros(1,total_reg);
 
 % Computing the Region length
 for pp=1: total_reg
     for dd=2:row-1
         for jj=2:col-1
             if(R(dd,jj)==pp)
             if((R(dd,jj)~=R(dd-1,jj))||(R(dd,jj)~=R(dd,jj-1))||(R(dd,jj)~=R(dd,jj+1))||(R(dd,jj)~=R(dd+1,jj+1)))
                 val=countregion(1,pp);
                 val=val+1;
                 countregion(1,pp)=val;
                 
             end
             end
         end
     end
     
 end
 
 
 %Computing the Shared boundary
 common_region=zeros(total_reg);
 for ss=2: total_reg-1
     for ff=2:row-1
         for jj=2:col-1
             if(R(ff,jj)==ss)
                if((R(ff,jj)~=R(ff-1,jj)))
                       regionNum=R(ff-1,jj);
                       if(regionNum~=0 && regionNum~=ss)
                                val=common_region(ss,regionNum);
                                val=val+1;
                               common_region(ss,regionNum)= val;
                       end
                end
                if((R(ff,jj)~=R(ff,jj-1)))
                       regionNum=R(ff,jj-1);
                            if(regionNum~=0 && regionNum~=ss)                           
                                val=common_region(ss,regionNum);
                                val=val+1;
                               common_region(ss,regionNum)= val;
                            end
                        
                end
                if((R(dd,jj)~=R(ff+1,jj+1)))
                        regionNum=R(ff+1,jj+1);
                        if(regionNum~=0 && regionNum~=ss) 
                            val=common_region(ss,regionNum);
                            val=val+1;
                            common_region(ss,regionNum)= val;
                       end
                end
                if((R(ff,jj)~=R(ff,jj+1)))
                       regionNum=R(ff,jj+1);
                        if(regionNum~=0 && regionNum~=ss)                         
                                val=common_region(ss,regionNum);
                                val=val+1;
                               common_region(ss,regionNum)= val;                      
                        end
                end
             end
         end
     end
     
 end
 

 
S=zeros(total_reg);
S_con=zeros(1,total_reg);
min_length=zeros(1,total_reg);
computedVal=zeros(1,total_reg);
k=256*256/total_reg;
   for i=1:total_reg-1
     %Computing Similarity score S_sim
     meandiff=abs(regionmeans(1,i)-regionmeans(1,i+1));
     array1=[1, regionstdeviation(1,i)+regionstdeviation(1,i+1)];
     max_std = max(array1);
     S_sim=meandiff/max_std;
     
     %Computing size constraint S_size
     array2 = [ regioncount(1,i), regioncount(1,i+1) ];
     mincacb=min(array2);
     array3 = [ 2, mincacb/k];
     S_size=min(array3);
     
     %Computing connectivity measure S_con
        S_con(1,i)=0.5;
      if(common_region(i,i+1)>0)
         minVal=[countregion(i), countregion(i+1)];
         min_length(1,i)=min(minVal);
         computedVal(1,i)=min_length(1,i)/(4*common_region(i,i+1));
         if((computedVal(1,i)>=0.5))
             if((computedVal(1,i)<=2.0))
                S_con(1,i)=computedVal(1,i);
             end
             elseif(computedVal(1,i)<0.5)
                     S_con(1,i)=0.5;
             else
                  S_con(1,i)= 2.0;
         end
     
        else
           S_con(1,i)=2.0;
       end
     %Computing final S(i,j) merge score between region i,j
     S(i,i+1)= ((S_sim ) * (sqrt(S_size)) )*  S_con(1,i);    
   end  
   %merge Score
   meanOfS=mean2(S);
   %region merging
   mergethreshold=meanOfS;      %Value of merge Score
    I2=I;
    intensity=zeros(1,total_reg);
   for i=1:total_reg-1
       if(S(i,i+1)>=0 && S(i,i+1)<= mergethreshold)
           %Merge Regions
            for ii=2:row-1
                for jj=2:col-1
                    if(R(ii,jj)==i)
                        intensity_val=intensity(1,i);
                        intensity_val=intensity_val+I(ii,jj);
                        intensity(1,i)=intensity_val;
                        
                    end
                
                end
            end
       end
   end
  
   
   for i=1:total_reg-1
       for i1=2:row-1
           for j1=2:col-1
               if(R(i1,j1)==i)
                   
                   I2(i1,j1)=intensity(1,i)/regioncount(1,i);
               end
           end
       end
   end
   
   figure(5),imshow(I2);
   title('Region Merging');
   
%-------------------Boundary Elimination--------------------------------%
I3=im2double(I2);
  global_thresh=2 ;
 
  a=1;
  b=1;
  c=1;
  d=1;
  e=2;
  f=2; 
  p=1;
  q=1;
  r=1;
  [row, col]=size(common_region);
   for i = 1:row
    for j=1:col
      if(common_region(i,j)>0)
          
           
        a=i;
        b=j;
      end
      [s,t]=size(regionNum);
         for k = 1:s
             for l=1:t
               if(regionNum(k,l) == i)
                   c=k;
                   d=l;
                   
                   e=I3(c,d);
                  
               end
               if(regionNum(k,l) == j)
                   p=k;
                   q=l;
                   
                   r=I3(p,k);
                  
                   
               end
             end
         end   
        
               if(e>r)
                 bondeli=((e-r)*common_region(i,j))/common_region(i,j);
                 if(bondeli < global_thresh)
                     I3(p,q)=I3(c,d);
                 
                 end
                     
               else 
                    bondeli=((r-e)*common_region(i,j))/common_region(i,j); 
                   
                     if(bondeli < global_thresh)
                       
                       I3(c,d)=I3(p,q);
                       
                    end
                  
               end
               
             
               
    end        
   end   

   figure(6);
   imshow(I3);
   title('Boundary Elimination');
   display('done');






 
 