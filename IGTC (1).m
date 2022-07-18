figure(1);
clf;

nodes =input("Enter no. of nodes: ");
k=input("Enter the value of k: ");

%% no grids
temp=0;% No.of grids temp x temp
if k==1
    temp=5;
elseif k==2
    temp=3;
elseif k==3
    temp=2;
end



%hold on;
field =input("Enter field size: ");
Rmax=100;

a=[];%x-coordinate of sensorNodes
b=[];%y-coordinate of sensorNodes
Wmatrix=zeros(nodes,nodes);

x_loc = rand(1,nodes)*field;
y_loc = rand(1,nodes)*field;

grid_size=k*(sqrt(2)*Rmax);

%% Deployment and connection of sensorNode to sensorNode
for i = 1:nodes
    
    a=[a,x_loc(i)];
    b=[b,y_loc(i)];

    final=[a; b];
    grid on
    set(gca,'XTick',0:grid_size:field);
    set(gca,'YTick',0:grid_size:field);
    plot(x_loc(i),y_loc(i),'*', 'Color','b', 'linewidth',2);
    text(x_loc(i)+1, y_loc(i)+1, num2str(i));
    hold on;
    %pause(0.5);
end 

%% Deployment of superNodes at MEAN position
no_of_superNodes=0;
sx=[];
sy=[];
lower_x=0;
lower_y=0;
coord_x=(grid_size);
coord_y=(grid_size);
for y=1:temp
  for x=1:temp 
    mean_x=[];
    mean_y=[];
   
    for j=1:nodes

        if (lower_x<=a(j)&&a(j)<=coord_x && lower_y<=b(j) && b(j)<=coord_y)
              mean_x=[mean_x,round(a(j))];
              mean_y=[mean_y,round(b(j))];
              M=[mean_x;mean_y];
             
                
        end

      
    end
     M=mean(M,2);
    %disp(M(1,1));
    
     
     if(~(M(1,1)==0 && M(2,1)==0))
        sx=[sx,M(1,1)];
        sy=[sy,M(2,1)];
       
        no_of_superNodes=no_of_superNodes+1;
     end
 
     
    lower_x=lower_x+grid_size;
    coord_x=coord_x+grid_size;
     M=zeros(2,1);
  end
  lower_x=0;
  lower_y=lower_y+grid_size;
  coord_x=grid_size;
  coord_y=coord_y+grid_size;
end
   

for i=1:no_of_superNodes
    
    plot(sx(i),sy(i),'d','Color','r', 'linewidth',4);
    hold on;
end
%% Connection of superNodes to sensorNodes
Smatrix=zeros(nodes,no_of_superNodes);
for i=1:nodes
    for j=1:no_of_superNodes
        dist1=sqrt((sx(j)-a(i))^2+(sy(j)-b(i))^2);
        if dist1<=Rmax
            Smatrix(i,j)=dist1;
            line([sx(j) a(i)], [sy(j) b(i)],'Color','r', 'LineStyle', '-');
             %matrix(j,j)=0;
        else
            Smatrix(i,j)=0;
        end
    end
    hold on;
end
disp("THis is Smatrix")

%% sensor to sensor node connection


for i=1:nodes
   
    for j = 1:nodes
        distance = sqrt((a(i) - a(j))^2 + (b(i) - b(j))^2);
        if distance <= Rmax
            Wmatrix(i, j) = distance; %there is a link
           
            line([a(i) a(j)], [b(i) b(j)],'Color','blue', 'LineStyle', ':');
            %Wmatrix(j,j)=0;
        else
            Wmatrix(i, j) = 0;
        end
    end
    hold on;
end



%% this Block of code is for K disjoint paths

S_transpose=transpose(Smatrix);
Wmatrix=horzcat(Wmatrix,Smatrix);

B=zeros(no_of_superNodes,no_of_superNodes);
for i=1:no_of_superNodes
    for j=1:no_of_superNodes
        if i==j
            B(i,j)=1;
        end
    end
end

S_transpose=horzcat(S_transpose,B);
final=vertcat(Wmatrix,S_transpose);

G=graph(final);
plot(G);

degree=zeros(nodes,no_of_superNodes); 
count=0;


%% this is the temporary data-structure to show the next hop of any sensor node


temporary=zeros(nodes,4);

% sensor nodes directly connected to supernodes
for i=1:nodes
    for j=1:no_of_superNodes
        temporary(i,1)=i;
        if (Smatrix(i,j)~=0)
            temporary(i,3)=sx(j);
            temporary(i,4)=sy(j);
            temporary(i,2)= j+nodes;
        end
    end     
end

% sensor nodes connected to supernodes through a max of 3 hops
for i=1:nodes
    for j=nodes+1:nodes+no_of_superNodes
        paths=allpaths(G,i,j,'MaxPathLength',3);
        if (~(isempty(paths)))
            if temporary(i,2)==0 && temporary(i,3)==0
                paths2=shortestpath(G,i,j);
                
                temporary(i,3)=a(paths2(2));
                temporary(i,4)=b(paths2(2));
                temporary(i,2)=paths2(2);
                
            end
        end
         
     end
 end

disp("no of supernodes: ");
disp(no_of_superNodes);

T = array2table(temporary,'VariableNames',{'Nodes','Next-Hop','X-coordinate','Y-coordinate'});
disp(T)


