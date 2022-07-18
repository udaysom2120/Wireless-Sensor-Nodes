nodes =input("Enter no. of nodes: ");
%firstfig =figure(1);
figure(1);
clf;
%subplot(1,2,1);
hold on;
field =input("Enter field size: ");
Rmax=100;
a=[];
b=[];

x_loc = rand(1,nodes)*field;
y_loc = rand(1,nodes)*field;
pause(3);

%% sensor node deployment and sensor node to sensor node connection
Wmatrix=zeros(nodes,nodes);

for i = 1:nodes
    plot(x_loc(i), y_loc(i), '*', 'Color','b', 'linewidth',2);
    text(x_loc(i)+1, y_loc(i)+1, num2str(i));
    a=[a,x_loc(i)];
    b=[b,y_loc(i)];
    final=[a; b];
    grid on
   
    set(gca,'XTick',0:100:field);
    set(gca,'YTick',0:100:field);
     for j = 1:nodes
        distance = sqrt((x_loc(i) - x_loc(j))^2 + (y_loc(i) - y_loc(j))^2);
        if distance <= Rmax
            Wmatrix(i, j) = 1; %there is a link
            line([x_loc(i) x_loc(j)], [y_loc(i) y_loc(j)],'Color','blue', 'LineStyle', ':');
            Wmatrix(j,j)=0;
        else
            Wmatrix(i, j) = 0;
        end
    end
    hold on;
   
end 


N = hist3(final','Ctrs',{10:100:field 10:100:field});
density_matrix = N.';


row=6;
column=6;
matrix9 = [0,0;0,0];
k=1;
for i=1:2:row
    for j=1:2:column
        B = density_matrix(i:i+1,j:j+1);
        matrix9(:,:,k) =B();
        k=k+1;
        %subMatrixCell{k}=subMatrix;
        % if you want to save all the
    end
end

mat=[0 200 400 0 200 400 0 200 400;
    0 0 0 200 200 200 400 400 400];
no_of_superNodes=0;

%figure(2);
%(1,2,2);
sx=[];
sy=[];
for i=1:9
    grid on;
    set(gca,'XTick',0:100:field);
    set(gca,'YTick',0:100:field);
    if nnz(matrix9(:,:,i)) == 4
        u=100+mat(1,i);
        v=100+mat(2,i);
        plot(u,v,'d', 'Color','k', 'linewidth',4);
        sx=[sx,u];
        sy=[sy,v];
        no_of_superNodes=no_of_superNodes+1;
        text(u+1,v+1,strcat('s',num2str(no_of_superNodes)))
        hold on;

    elseif nnz(matrix9(:,:,i)) == 3
        u=100+mat(1,i);
        v=100+mat(2,i);
        
        plot(u,v,'d', 'Color','k', 'linewidth',4);
        sx=[sx,u];
        sy=[sy,v];
        no_of_superNodes=no_of_superNodes+1;
        text(u+1,v+1,strcat('s',num2str(no_of_superNodes)))
        hold on;

    elseif nnz(matrix9(:,:,i)) == 2
        if matrix9(1,1,i) == 0
            if matrix9(1,2,i) == 0
                u=100+mat(1,i);
                v=150+mat(2,i);
               
                plot(u,v,'d', 'Color','k', 'linewidth',4);
                sx=[sx,u];
                sy=[sy,v];
                no_of_superNodes=no_of_superNodes+1;
                text(u+1,v+1,strcat('s',num2str(no_of_superNodes)))
                hold on;

            elseif matrix9(2,1,i) == 0
                u=150+mat(1,i);
                v=100+mat(2,i);
              
                plot(u,v,'d', 'Color','k', 'linewidth',4);
                sx=[sx,u];
                sy=[sy,v];
                no_of_superNodes=no_of_superNodes+1;
                text(u+1,v+1,strcat('s',num2str(no_of_superNodes)))
                hold on;

            elseif matrix9(2,2,i) == 0
                u=100+mat(1,i);
                v=100+mat(2,i);
                plot(u,v,'d', 'Color','k', 'linewidth',4);
                sx=[sx,u];
                sy=[sy,v];
                no_of_superNodes=no_of_superNodes+1;
                text(u+1,v+1,strcat('s',num2str(no_of_superNodes)))
                hold on;
            end

        elseif matrix9(2,1,i) == 0
            if matrix9(1,1,i) == 0
                u=150+mat(1,i);
                v=100+mat(2,i);
               
                plot(u,v,'d', 'Color','k', 'linewidth',4);
                sx=[sx,u];
                sy=[sy,v];
                hold on;
                no_of_superNodes=no_of_superNodes+1;
                text(u+1,v+1,strcat('s',num2str(no_of_superNodes)));

            elseif matrix9(1,2,i) == 0
                u=100+mat(1,i);
                v=100+mat(2,i);
                %disp("C3.2.2");
                %disp(i);
                %disp(u);disp(v);
                plot(u,v,'d', 'Color','k', 'linewidth',4);
                sx=[sx,u];
                sy=[sy,v];
                hold on;
                no_of_superNodes=no_of_superNodes+1;
                text(u+1,v+1,strcat('s',num2str(no_of_superNodes)));

            elseif matrix9(2,2,i) == 0
                u=100+mat(1,i);
                v=50+mat(2,i);
               
                plot(u,v,'d', 'Color','k', 'linewidth',4);
                sx=[sx,u];
                sy=[sy,v];
                no_of_superNodes=no_of_superNodes+1;
                text(u+1,v+1,strcat('s',num2str(no_of_superNodes)));
                hold on;
            end

        elseif matrix9(2,2,i) == 0
            if matrix9(1,1,i) == 0
                u=100+mat(1,i);
                v=100+ mat(2,i);
              
                plot(u,v,'d', 'Color','k', 'linewidth',4);
                sx=[sx,u];
                sy=[sy,v];
                no_of_superNodes=no_of_superNodes+1;
                text(u+1,v+1,strcat('s',num2str(no_of_superNodes)));
                hold on;

            elseif matrix9(1,2,i) == 0
                u=50+mat(1,i);
                v=100+ mat(2,i);
              
                plot(u,v,'d', 'Color','k', 'linewidth',4);
                sx=[sx,u];
                sy=[sy,v];
                no_of_superNodes=no_of_superNodes+1;
                text(u+1,v+1,strcat('s',num2str(no_of_superNodes)));
                hold on;

            elseif matrix9(2,1,i) == 0
                u=100+mat(1,i);
                v=50+ mat(2,i);
               
                plot(u,v,'d', 'Color','k', 'linewidth',4);
                sx=[sx,u];
                sy=[sy,v];
                no_of_superNodes=no_of_superNodes+1;
                text(u+1,v+1,strcat('s',num2str(no_of_superNodes)));
                hold on;
            end
        end

    elseif nnz(matrix9(:,:,i)) == 1
        if matrix9(1,1,i) ~= 0
            u=50+mat(1,i);
            v=50+mat(2,i);
         
            plot(u,v,'d', 'Color','k', 'linewidth',4)
            sx=[sx,u];
            sy=[sy,v];
            no_of_superNodes=no_of_superNodes+1;
            text(u+1,v+1,strcat('s',num2str(no_of_superNodes)));
            hold on;

        elseif matrix9(1,2,i) ~=0
            u=150+mat(1,i);
            v=50+mat(2,i);
         
            plot(u,v,'d', 'Color','k', 'linewidth',4);
            sx=[sx,u];
            sy=[sy,v];
            no_of_superNodes=no_of_superNodes+1;
            text(u+1,v+1,strcat('s',num2str(no_of_superNodes)));
            hold on;

        elseif matrix9(2,1,i)~=0
            u=50+mat(1,i);
            v=150+mat(2,i);
          
            plot(u,v,'d', 'Color','k', 'linewidth',4);
            sx=[sx,u];
            sy=[sy,v];
            no_of_superNodes=no_of_superNodes+1;
            text(u+1,v+1,strcat('s',num2str(no_of_superNodes)));
            hold on;

        elseif matrix9(2,2,i)~=0
            u=150+mat(1,i);
            v=150+mat(2,i);
         
            plot(u,v,'d', 'Color','k', 'linewidth',4);
            sx=[sx,u];
            sy=[sy,v];
            no_of_superNodes=no_of_superNodes+1;
            text(u+1,v+1,strcat('s',num2str(no_of_superNodes)));
            hold on;
        end
    end
end

%% sensor node to supernode connection
Smatrix=zeros(nodes,no_of_superNodes);
Smax=100;
%CM = jet(nodes);  % See the help for COLORMAP to see other choices.
for i=1:nodes
    for j=1:no_of_superNodes
        dist=sqrt((sx(j)-a(i))^2+(sy(j)-b(i))^2);
        if dist<=Smax
            Smatrix(i,j)=1;
            line([sx(j) a(i)], [sy(j) b(i)],'Color','r', 'LineStyle', ':');
             %matrix(j,j)=0;
        else
            Smatrix(i,j)=0;
        end
    end
end


%% for k disjoint paths
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
final2=vertcat(Wmatrix,S_transpose);

G=graph(final2);


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



