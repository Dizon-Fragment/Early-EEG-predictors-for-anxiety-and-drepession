classdef bubbleDigraph_dgz
% @author : Dizonz
% 公众号  : slandarer随笔 
% Zhaoxu Liu / Dizonz (2023). bubble digraph 
% (https://www.mathworks.com/matlabcentral/fileexchange/125140-bubble-digraph), 
% MATLAB Central File Exchange. 检索来源 2023/2/21.
    properties
        ax,arginList={'ColorOrder','Colormap','ClassName','NodeName','BubbleSize','RClass','RNode', 'Pthresh','LineWidth'}
        ClassName   % 类名称
        NodeName    % 节点名称
        % 节点配色
        ColorOrder=[0.6510    0.8078    0.8902;    0.1216    0.4706    0.7059;    0.6980    0.8745    0.5412
                    0.2000    0.6275    0.1725;    0.9843    0.6039    0.6000;    0.8902    0.1020    0.1098
                    0.9922    0.7490    0.4353;    1.0000    0.4980         0;    0.7922    0.6980    0.8392
                    0.4157    0.2392    0.6039;    1.0000    1.0000    0.6000;    0.6941    0.3490    0.1569];
        % 连线配色
        Colormap=[1.0000    0.9686    0.9529;    0.9980    0.9454    0.9307;    0.9960    0.9221    0.9084;    0.9939    0.8988    0.8861
                  0.9920    0.8750    0.8630;    0.9910    0.8477    0.8336;    0.9900    0.8204    0.8043;    0.9890    0.7930    0.7750;    
                  0.9877    0.7629    0.7502;    0.9857    0.7245    0.7390;    0.9837    0.6860    0.7279;    0.9817    0.6476    0.7168;    
                  0.9793    0.6027    0.7022;    0.9762    0.5470    0.6820;    0.9732    0.4913    0.6617;    0.9701    0.4357    0.6415;    
                  0.9555    0.3815    0.6263;    0.9292    0.3289    0.6162;    0.9028    0.2763    0.6061;    0.8765    0.2237    0.5960;    
                  0.8369    0.1717    0.5763;    0.7894    0.1201    0.5510;    0.7418    0.0684    0.5257;    0.6942    0.0168    0.5004;    
                  0.6429    0.0039    0.4888;    0.5903    0.0039    0.4817;    0.5376    0.0039    0.4746;    0.4850    0.0039    0.4676;    
                  0.4350    0.0030    0.4552;    0.3855    0.0020    0.4420;    0.3359    0.0010    0.4288;    0.2863         0    0.4157];
        Data,Class          % 有向图及节点分类
        BubbleSize=[5,25];  % 节点气泡大小范围[min,max]
        LineWidth=1.5;      % 连接线粗细
        AlphaLim=[.1,.9];
        ClassSet,ClassNum,
        RClass=1.25;
        RNode=1.08;
        ColorList;
        Pthresh
        bubbleHdl,nodeLabelHdl,classLabelHdl
    end

    methods
        function obj=bubbleDigraph_dgz(Data,Class,varargin)
            obj.Data=Data;
            obj.Class=Class(:);
            obj.ClassSet=unique(Class);
            obj.ClassNum=length(obj.ClassSet);

            for i=1:size(obj.Data,1)
                obj.NodeName{i}='';
            end
            for i=1:obj.ClassNum
                obj.ClassName{i}='';
            end
            
            % 获取其他数据
            disp(char([64 97 117 116 104 111 114 32 58 32,...
                 115 108 97 110 100 97 114 101 114]))
            for i=1:2:(length(varargin)-1)
                tid=ismember(obj.arginList,varargin{i});
                if any(tid)
                    obj.(obj.arginList{tid})=varargin{i+1};
                end
            end
            if obj.ClassNum>size(obj.ColorOrder,1)
                obj.ColorOrder=[obj.ColorOrder;rand([obj.ClassNum,3])];
            end
        end
        function obj=draw(obj)
            obj.ax=gca;hold on;
            colormap(obj.Colormap)
            obj.ax.XLim=[-1.2,1.2];
            obj.ax.YLim=[-1.2,1.2];
            obj.ax.XTick=[];
            obj.ax.YTick=[];
            obj.ax.XColor='none';
            obj.ax.YColor='none';
            obj.ax.PlotBoxAspectRatio=[1,1,1];

            % 调整初始界面大小
            fig=obj.ax.Parent;
            fig.Color=[1,1,1];
            if max(fig.Position(3:4))<600
                fig.Position(3:4)=1.8.*fig.Position(3:4);
                fig.Position(1:2)=fig.Position(1:2)./3;
            end
            % draw bubble
            thetaList=linspace(0,2*pi,size(obj.Data,1)+1);thetaList(end)=[];
            XList=cos(thetaList);YList=sin(thetaList);
            obj.bubbleHdl=bubblechart(XList,YList,sum(abs(obj.Data)));
            bubblesize(obj.BubbleSize)

            obj.ColorList=zeros(size(obj.Data,1),3);
            tmp_idx = cellfun(@isempty, obj.ClassName);
            for i=1:length(obj.ClassSet)
                if tmp_idx(i)
                    obj.ColorList(obj.Class==obj.ClassSet(i),:)=...
                    repmat([1,1,1],sum(obj.Class==obj.ClassSet(i)),1);
                else
                    obj.ColorList(obj.Class==obj.ClassSet(i),:)=...
                        repmat(obj.ColorOrder(i,:),sum(obj.Class==obj.ClassSet(i)),1);
                end
            end
            obj.bubbleHdl.CData=obj.ColorList;
            % 绘制连线
            alphaData=abs(obj.Data);
            alphaData=alphaData-min(min(alphaData));
            alphaData=alphaData./max(max(alphaData));
            alphaData=alphaData.*(obj.AlphaLim(2)-obj.AlphaLim(1))+obj.AlphaLim(1);
            for i=1:size(obj.Data,1)
                for j=1:i-1
                    if obj.Data(i,j)~=0
                        bezierX=[cos(thetaList(i)),0,cos(thetaList(j))].*.93;
                        bezierY=[sin(thetaList(i)),0,sin(thetaList(j))].*.93;
                        bezierPnts=bezierCurve([bezierX',bezierY'],100);
                        bezierX=[bezierPnts(:,1);nan];
                        bezierY=[bezierPnts(:,2);nan];
                        if obj.Data(i,j) < obj.Pthresh
                            fill(bezierX,bezierY,linspace(-1,1,101).*obj.Data(i,j)./abs(obj.Data(i,j)),'EdgeColor',[0.5 0.5 0.5],...
                            'LineWidth',obj.LineWidth,'EdgeAlpha',0.2)
                        else
                            fill(bezierX,bezierY,linspace(-1,1,101).*obj.Data(i,j)./abs(obj.Data(i,j)),'EdgeColor','interp',...
                                'LineWidth',obj.LineWidth,'EdgeAlpha',alphaData(i,j))
                        end
                    end
                end
            end
            % 绘制节点标签
            for i=1:size(obj.Data,1)
                Ti=thetaList(i);
                rotation=Ti/pi*180;
                if rotation>90&&rotation<270
                    rotation=rotation+180;
                    obj.nodeLabelHdl(i)=text(cos(Ti).*obj.RNode,sin(Ti).*obj.RNode,obj.NodeName{i},...
                        'Rotation',rotation,'HorizontalAlignment','right','FontSize',9);
                else
                    obj.nodeLabelHdl(i)=text(cos(Ti).*obj.RNode,sin(Ti).*obj.RNode,obj.NodeName{i},...
                        'Rotation',rotation,'FontSize',9);
                end
            end
            % 绘制类标签
            for i=1:obj.ClassNum
                CTi=mean(thetaList(obj.Class==obj.ClassSet(i)));
                rotation=CTi/pi*180;
                if rotation>0&&rotation<180
                    obj.classLabelHdl(i)=text(cos(CTi).*obj.RClass,sin(CTi).*obj.RClass,obj.ClassName{i},'FontSize',14,'FontName','Arial',...
                    'HorizontalAlignment','center','Rotation',-(.5*pi-CTi)./pi.*180);
                else
                    obj.classLabelHdl(i)=text(cos(CTi).*obj.RClass,sin(CTi).*obj.RClass,obj.ClassName{i},'FontSize',14,...
                    'HorizontalAlignment','center','Rotation',-(1.5*pi-CTi)./pi.*180);
                end
            end
            % 贝塞尔函数
            function pnts=bezierCurve(pnts,N)
                t=linspace(0,1,N);
                p=size(pnts,1)-1;
                coe1=factorial(p)./factorial(0:p)./factorial(p:-1:0);
                coe2=((t).^((0:p)')).*((1-t).^((p:-1:0)'));
                pnts=(pnts'*(coe1'.*coe2))';
            end
        end
        % 修改气泡颜色
        function obj=setBubbleColor(obj,ColorList)
            obj.ColorOrder=ColorList;
            tmp_idx = cellfun(@isempty, obj.ClassName);
            for i=1:length(obj.ClassSet)
                if tmp_idx(i)
                    obj.ColorList(obj.Class==obj.ClassSet(i),:)=...
                    repmat([1,1,1],sum(obj.Class==obj.ClassSet(i)),1);
                else
                    obj.ColorList(obj.Class==obj.ClassSet(i),:)=...
                        repmat(obj.ColorOrder(i,:),sum(obj.Class==obj.ClassSet(i)),1);
                end
            end
            set(obj.bubbleHdl,'CData',obj.ColorList);
        end
        % 修改气泡其他属性
        function obj=setBubble(obj,varargin)
            set(obj.bubbleHdl,varargin{:});
        end
        % 设置标签
        function setNodeLabel(obj,varargin)
            for i=1:size(obj.Data,1)
                set(obj.nodeLabelHdl(i),varargin{:});
            end
        end

        function setClassLabel(obj,varargin)
            for i=1:obj.ClassNum
                set(obj.classLabelHdl(i),varargin{:});
            end
        end
    end
% @author : Dizonz
% 公众号  : Dizonz随笔 
% Zhaoxu Liu / slandarer (2023). bubble digraph 
% (https://www.mathworks.com/matlabcentral/fileexchange/125140-bubble-digraph), 
% MATLAB Central File Exchange. 检索来源 2023/2/21.
end