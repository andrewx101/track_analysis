classdef GSER_dlg_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        GSERUIFigure        matlab.ui.Figure
        bpgSmoothAlgorithm  matlab.ui.container.ButtonGroup
        optSLM              matlab.ui.control.RadioButton
        optRLowess          matlab.ui.control.RadioButton
        optMovingAvg        matlab.ui.control.RadioButton
        btnOK               matlab.ui.control.Button
        sldSpan             matlab.ui.control.Slider
        SpanLabel           matlab.ui.control.Label
        sldMaxDt            matlab.ui.control.Slider
        MaxLabel            matlab.ui.control.Label
        sldMinDt            matlab.ui.control.Slider
        MinLabel            matlab.ui.control.Label
        axePlot2            matlab.ui.control.UIAxes
        axePlot1            matlab.ui.control.UIAxes
    end


    properties (Access = private)
        callingApp
        msd
        T
        a

        min_dt
        max_dt

        msd_cut
        msd_smooth
        span
        G
        smoothAlgorithm=1
    end

    methods (Access = private)

        function results = recal(app)




            % For modes 2 & 3 only: interpolate into even spacing under log scale
            N=length(app.msd(:,1));
            dt=logspace(log10(app.msd(1,1)),log10(app.msd(end,1)),N);
            dt=dt(:);
            msd_raw=[dt,interp1(app.msd(:,1),app.msd(:,2),dt,'linear','extrap'),...
                interp1(app.msd(:,1),app.msd(:,3),dt,'linear','extrap'),...
                interp1(app.msd(:,1),app.msd(:,4),dt,'linear','extrap')];

            dt=app.msd(:,1); % extract the raw dt of msd

            % smooth
            switch app.smoothAlgorithm
                case 1 %  moving-average with fixed logarithmic time-span
                    w=app.sldSpan.Value; % portion of a time decade
                    % In this mode, sldSpan.Value should lie within [0,1];

                    % pre-smooth the msd by moving average using this span
                    msd_smoothed=zeros(length(dt),1);
                    for j=1:length(dt)
                        dt1=dt(j).*10.^(-w./2);
                        dt2=dt(j).*10.^(w./2);
                        ind=find(dt>=dt1 & dt<=dt2);
                        if isempty(ind)
                            msd_smoothed(j)=app.msd(j,2);
                        else
                            msd_temp=10.^smooth(log10(app.msd(:,2)),length(ind),'moving');
                            msd_smoothed(j)=msd_temp(j);
                        end
                    end
                    app.msd_smooth=[dt(:),msd_smoothed(:)];

                    % Calculate the local slope alpha
                    alpha=gradient(log10(app.msd_smooth(:,2)))./gradient(log10(app.msd_smooth(:,1)));

                    % smooth alpha again by the same method
                    alpha_smoothed=zeros(length(dt),1);
                    for j=1:length(dt)
                        dt1=dt(j).*10.^(-w./2);
                        dt2=dt(j).*10.^(w./2);
                        ind=find(dt>=dt1 & dt<=dt2);
                        if isempty(ind)
                            alpha_smoothed(j)=alpha(j);
                        else
                            alpha_temp=smooth(alpha,length(ind),'moving');
                            alpha_smoothed(j)=alpha_temp(j);
                        end
                    end
                    alpha=alpha_smoothed;

                case 2
                    app.msd_smooth=[msd_raw(:,1),10.^...
                        smooth(log10(msd_raw(:,1)),log10(msd_raw(:,2)),...
                        app.sldSpan.Value,'rloess')];
                    alpha=gradient(log10(app.msd_smooth(:,2)))./gradient(log10(app.msd_smooth(:,1)));
                    alpha=smooth(app.msd(:,1),alpha,app.sldSpan.Value,'rloess');
                    dt=msd_raw(:,1);
                case 3
                    slm=slmengine(log10(msd_raw(:,1)),log10(msd_raw(:,2)),'increasing','on',...
                        'maxslope',1,'minslope',0,'C2','on','knots',app.sldSpan.Value);
                    app.msd_smooth=[msd_raw(:,1),10.^slmeval(log10(msd_raw(:,1)),slm)];
                    alpha=gradient(log10(app.msd_smooth(:,2)))./gradient(log10(app.msd_smooth(:,1)));
                    slm=slmengine(log10(app.msd(:,1)),alpha,'maxvalue',0.999363380439839,...
                        'minvalue',0.000636619560161118,'knots',app.sldSpan.Value,'C2','on');
                    alpha=slmeval(log10(app.msd(:,1)),slm);
                    dt=msd_raw(:,1);
            end

            % cut

            [~,ind1]=min(abs(dt(:,1)-10^app.sldMinDt.Value));
            [~,ind2]=min(abs(dt(:,1)-10^app.sldMaxDt.Value));

            app.msd_cut=app.msd_smooth(ind1:ind2,:);
            alpha=alpha(ind1:ind2);

            % GSER

            w=1./app.msd_cut(1:end,1);
            Ga=1.38e-23*(app.T+273.15)./(pi*app.a.*app.msd_cut(1:end,2).*gamma(1+alpha));
            Gp=abs(Ga).*cos(pi.*alpha./2);
            Gpp=abs(Ga).*sin(pi.*alpha./2);

            G=[w,Gp,Gpp];
            G(G(:,2)<0 | G(:,3)<0 | isnan(G(:,2)) | isnan(G(:,3)),:)=[];
            app.G=G;

            % plot MSD
            cla(app.axePlot1,'reset')
            reset(app.axePlot1)
            app.axePlot1.Box='on';
            loglog(app.axePlot1,app.msd(:,1),app.msd(:,2),'.');

            app.axePlot1.XLabel.String='lag time (s)';
            app.axePlot1.YLabel.String='MSD (m^2)';
            app.axePlot1.XLim=10.^[floor(log10(min(app.msd(:,1)))),ceil(log10(max(app.msd(:,1))))];
            app.axePlot1.YLim=10.^[floor(log10(min(app.msd(:,2)))),ceil(log10(max(app.msd(:,2))))];
            app.axePlot1.XTick=10.^(floor(log10(min(app.msd(:,1)))):ceil(log10(max(app.msd(:,1)))));
            app.axePlot1.YTick=10.^(floor(log10(min(app.msd(:,2)))):ceil(log10(max(app.msd(:,2)))));

            hold(app.axePlot1)
            plot(app.axePlot1,app.msd_cut(:,1),app.msd_cut(:,2),'LineWidth',1)
            plot(app.axePlot1,[10^app.sldMinDt.Value,10^app.sldMinDt.Value],app.axePlot1.YLim,'k--');
            plot(app.axePlot1,[10^app.sldMaxDt.Value,10^app.sldMaxDt.Value],app.axePlot1.YLim,'k--');
            hold(app.axePlot1)

            % plot GSER
            cla(app.axePlot2,'reset');

            if ~isempty(app.G)

                app.axePlot2.Box='on';
                loglog(app.axePlot2,app.G(:,1),app.G(:,2:3),'LineWidth',1);

                app.axePlot2.XLabel.String='angular frequency (rad/s)';
                app.axePlot2.YLabel.String='moduli (Pa)';
                app.axePlot2.XLim=10.^[floor(log10(min(app.G(:,1)))),ceil(log10(max(app.G(:,1))))];

                %             app.axePlot2.YLim=10.^[floor(log10(min((app.G(:,2:3)),[],'all'))),ceil(log10(max(app.G(:,2:3),[],'all')))];
                %             app.axePlot2.YLim=1.38e-23.*(app.T+273.15)./...
                %                 [app.axePlot1.YLim(2).*app.a.*app.axePlot2.XLim(2).*pi,...
                %                 app.axePlot1.YLim(1).*app.a.*app.axePlot2.XLim(1).*pi];



                app.axePlot2.XTick=10.^(floor(log10(min(app.G(:,1)))):ceil(log10(max(app.G(:,1)))));
                axis(app.axePlot2,'auto y');
                lim=axis(app.axePlot2);
                app.axePlot2.YLim=10.^[floor(log10(lim(3))),ceil(log10(lim(4)))];
                min_y=1e-308;
                
                if app.axePlot2.YLim(1)==0
                    app.axePlot2.YLim(1)=min_y;
                end
                if isinf(app.axePlot2.YLim(2))
                    app.axePlot2.YLim(2)=realmax;
                end
                % axis(app.axePlot2,'auto y');
                yticks=10.^(floor(log10(app.axePlot2.YLim(1))):1:ceil(log10(app.axePlot2.YLim(2))));
                s=1;
                while length(yticks)>7
                    yticks=10.^(floor(log10(app.axePlot2.YLim(1))):s:ceil(log10(app.axePlot2.YLim(2))));
                    s=s+1;
                end

                app.axePlot2.YTick=yticks;


                lg=legend(app.axePlot2);

                lg.Interpreter='latex';
                lg.String=[{'$G\prime$'},{'$G\prime\prime$'}];
                lg.Location='best';
                lg.Box='off';
            else
                tx=text(app.axePlot2,app.axePlot2.XLim(1),...
                    app.axePlot2.YLim(1)+(app.axePlot2.YLim(2)-app.axePlot2.YLim(1))/2,...
                    {'No valid results.','Try smaller span.'});
                tx.VerticalAlignment='middle';

            end



        end
    end


    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app, caller, msd, T, a)
            app.callingApp=caller;
            app.msd=msd;
            app.T=T;
            app.a=a;
            app.min_dt=app.msd(1,1);
            app.max_dt=app.msd(end,1);

            % update sliders
            app.sldMinDt.Limits=log10([min(app.msd(:,1)),max(app.msd(:,1))]);
            app.sldMaxDt.Limits=app.sldMinDt.Limits;
            decade_points=ceil(log10(min(app.msd(:,1)))):floor(log10(max(app.msd(:,1))));
            app.sldMinDt.MajorTicks=decade_points;
            app.sldMaxDt.MajorTicks=decade_points;
            app.sldMinDt.MajorTickLabels=arrayfun(@(i) num2str(10^decade_points(i)),1:length(decade_points),'UniformOutput',0);
            app.sldMaxDt.MajorTickLabels=app.sldMinDt.MajorTickLabels;
            app.sldMinDt.MinorTicks=floor(log10(min(app.msd(:,1)))/0.2)*0.2:0.2:ceil(log10(max(app.msd(:,1)))/0.2)*0.2;
            app.sldMaxDt.MinorTicks=app.sldMinDt.MinorTicks;
            app.sldMinDt.Value=log10(min(app.msd(:,1)));
            app.sldMaxDt.Value=log10(max(app.msd(:,1)));

            app.sldSpan.Limits=[0,1];
            app.sldSpan.Value=0.2;
            drawnow

            app.recal;


            %             % pre-smooth
            %             app.recal
            %
            %
            %
            %
            %             % plot MSD
            %             reset(app.axePlot1)
            %             app.axePlot1.Box='on';
            %             loglog(app.axePlot1,app.msd(:,1),app.msd(:,2),'.');
            %
            %             app.axePlot1.XLabel.String='lag time (s)';
            %             app.axePlot1.YLabel.String='MSD (m^2)';
            %             app.axePlot1.XLim=10.^[floor(log10(min(app.msd(:,1)))),ceil(log10(max(app.msd(:,1))))];
            %             app.axePlot1.YLim=10.^[floor(log10(min(app.msd(:,2)))),ceil(log10(max(app.msd(:,2))))];
            %             app.axePlot1.XTick=10.^(floor(log10(min(app.msd(:,1)))):ceil(log10(max(app.msd(:,1)))));
            %             app.axePlot1.YTick=10.^(floor(log10(min(app.msd(:,2)))):ceil(log10(max(app.msd(:,2)))));
            %
            %             hold(app.axePlot1)
            %             plot(app.axePlot1,app.msd_cut(:,1),app.msd_cut(:,2),'LineWidth',1)
            %             plot(app.axePlot1,[app.min_dt,app.min_dt],app.axePlot1.YLim,'k--');
            %             plot(app.axePlot1,[app.max_dt,app.max_dt],app.axePlot1.YLim,'k--');
            %             hold(app.axePlot1)
            %
            %             % pre-process
            %
            %             msd_input=app.msd_cut;
            %
            %             alpha=gradient(log10(msd_input(:,2)))./gradient(log10(msd(:,1)));
            %             w=1./msd_input(1:end,1);
            %
            %             Ga=1.38e-23*(app.T+273.15)./(pi*a.*msd_input(1:end,2).*gamma(1+alpha));
            %             Gp=Ga.*cos(pi.*alpha./2);
            %             Gpp=Ga.*sin(pi.*alpha./2);
            %
            %             app.G=[w,Gp,Gpp];
            %
            %
            %             % plot GSER
            %             cla(app.axePlot2,'reset');
            %
            %             app.axePlot2.Box='on';
            %             loglog(app.axePlot2,app.G(:,1),app.G(:,2:3),'LineWidth',1);
            %
            %             app.axePlot2.XLabel.String='angular frequency (rad/s)';
            %             app.axePlot2.YLabel.String='Moduli (Pa)';
            %             app.axePlot2.XLim=10.^[floor(log10(min(app.G(:,1)))),ceil(log10(max(app.G(:,1))))];
            %             app.axePlot2.YLim=1.38e-23.*(app.T+273.15)./...
            %                 [app.axePlot1.YLim(2).*app.a.*app.axePlot2.XLim(2).*pi,...
            %                 app.axePlot1.YLim(1).*app.a.*app.axePlot2.XLim(1).*pi];
            %
            %
            %             app.axePlot2.XTick=10.^(floor(log10(min(app.G(:,1)))):ceil(log10(max(app.G(:,1)))));
            %             app.axePlot2.YTick=10.^(floor(log10(app.axePlot2.YLim(1))):ceil(log10(app.axePlot2.YLim(2))));
            %
            %             app.axePlot2.YTick=10.^(floor(log10(min([app.G(app.G(:,2)>0,2);app.G(app.G(:,3)>0,3)]))):...
            %                 ceil(log10(max([app.G(app.G(:,2)>0,2);app.G(app.G(:,3)>0,3)]))));
            %             lg=legend(app.axePlot2);
            %
            %             lg.Interpreter='latex';
            %             lg.String=[{'$G\prime$'},{'$G\prime\prime$'}];
            %             lg.Location='best';






        end

        % Value changed function: sldMinDt
        function sldMinDtValueChanged(app, event)
            value = app.sldMinDt.Value;
            slider_value=10^value;
            [~,ind1]=min(abs(app.msd(:,1)-slider_value));
            app.sldMinDt.Value=log10(app.msd(ind1,1));

            [~,ind2]=min(abs(app.msd(:,1)-10^app.sldMaxDt.Value));
            if ind2-ind1<10 && ind1+10<length(app.msd(:,1))
                ind2=ind1+10;
                app.sldMaxDt.Value=log10(app.msd(ind2,1));
            elseif ind2-ind1<10 && ind1+10>=length(app.msd(:,1))
                app.sldMinDt.Value=log10(app.msd(end-10,1));
                app.sldMaxDt.Value=log10(app.msd(end,1));
            end
            drawnow

            %             app.msd_cut=app.msd_smooth(ind1:ind2,:);
            %             app.min_dt=app.msd(ind1,1);
            %             app.max_dt=app.msd(ind2,1);
            %
            %             reset(app.axePlot1)
            %             app.axePlot1.Box='on';
            %             loglog(app.axePlot1,app.msd(:,1),app.msd(:,2),'.');
            %
            %             app.axePlot1.XLabel.String='lag time (s)';
            %             app.axePlot1.YLabel.String='MSD (m^2)';
            %             app.axePlot1.XLim=10.^[floor(log10(min(app.msd(:,1)))),ceil(log10(max(app.msd(:,1))))];
            %             app.axePlot1.YLim=10.^[floor(log10(min(app.msd(:,2)))),ceil(log10(max(app.msd(:,2))))];
            %             app.axePlot1.XTick=10.^(floor(log10(min(app.msd(:,1)))):ceil(log10(max(app.msd(:,1)))));
            %             app.axePlot1.YTick=10.^(floor(log10(min(app.msd(:,2)))):ceil(log10(max(app.msd(:,2)))));
            %
            %             hold(app.axePlot1)
            %             plot(app.axePlot1,app.msd_cut(:,1),app.msd_cut(:,2),'LineWidth',1)
            %             plot(app.axePlot1,[app.min_dt,app.min_dt],app.axePlot1.YLim,'k--');
            %             plot(app.axePlot1,[app.max_dt,app.max_dt],app.axePlot1.YLim,'k--');
            %             hold(app.axePlot1)

            app.recal;


        end

        % Value changing function: sldMinDt
        function sldMinDtValueChanging(app, event)
            % changingValue = event.Value;
            % slider_value=10^changingValue;
            % [~,ind1]=min(abs(app.msd(:,1)-slider_value));
            %
            %
            % [~,ind2]=min(abs(app.msd(:,1)-10^app.sldMaxDt.Value));
            % if ind2-ind1<10 && ind1+10<length(app.msd(:,1))
            %     ind2=ind1+10;
            %     app.sldMaxDt.Value=log10(app.msd(ind2,1));
            % elseif ind2-ind1<10 && ind1+10>=length(app.msd(:,1))
            %     app.sldMinDt.Value=log10(app.msd(end-10,1));
            %     app.sldMaxDt.Value=log10(app.msd(end,1));
            % end
            % drawnow

        end

        % Value changed function: sldMaxDt
        function sldMaxDtValueChanged(app, event)
            value = app.sldMaxDt.Value;
            slider_value=10^value;


            [~,ind2]=min(abs(app.msd(:,1)-slider_value));

            app.sldMaxDt.Value=log10(app.msd(ind2,1));

            [~,ind1]=min(abs(app.msd(:,1)-10^app.sldMinDt.Value));
            if ind2-ind1<10&&ind2-10>=1
                ind1=ind2-10;
                app.sldMinDt.Value=log10(app.msd(ind1,1));
            elseif ind2-ind1<10&&ind2-10<1
                app.sldMinDt.Value=log10(app.msd(1,1));
                app.sldMaxDt.Value=log10(app.msd(11,1));
            end
            drawnow

            %             app.msd_cut=app.msd_smooth(ind1:ind2,:);
            %             app.min_dt=app.msd(ind1,1);
            %             app.max_dt=app.msd(ind2,1);
            %
            %
            %             reset(app.axePlot1)
            %             app.axePlot1.Box='on';
            %             loglog(app.axePlot1,app.msd(:,1),app.msd(:,2),'.');
            %
            %             app.axePlot1.XLabel.String='lag time (s)';
            %             app.axePlot1.YLabel.String='MSD (m^2)';
            %             app.axePlot1.XLim=10.^[floor(log10(min(app.msd(:,1)))),ceil(log10(max(app.msd(:,1))))];
            %             app.axePlot1.YLim=10.^[floor(log10(min(app.msd(:,2)))),ceil(log10(max(app.msd(:,2))))];
            %             app.axePlot1.XTick=10.^(floor(log10(min(app.msd(:,1)))):ceil(log10(max(app.msd(:,1)))));
            %             app.axePlot1.YTick=10.^(floor(log10(min(app.msd(:,2)))):ceil(log10(max(app.msd(:,2)))));
            %
            %             hold(app.axePlot1)
            %             plot(app.axePlot1,app.msd_cut(:,1),app.msd_cut(:,2),'LineWidth',1)
            %             plot(app.axePlot1,[app.min_dt,app.min_dt],app.axePlot1.YLim,'k--');
            %             plot(app.axePlot1,[app.max_dt,app.max_dt],app.axePlot1.YLim,'k--');
            %             hold(app.axePlot1)
            app.recal;

        end

        % Value changing function: sldMaxDt
        function sldMaxDtValueChanging(app, event)
            changingValue = event.Value;

            slider_value=10^changingValue;
            [~,ind2]=min(abs(app.msd(:,1)-slider_value));
            app.sldMaxDt.Value=log10(app.msd(ind2,1));

            [~,ind1]=min(abs(app.msd(:,1)-10^app.sldMinDt.Value));
            if ind2-ind1<10&&ind2-10>=1
                ind1=ind2-10;
                app.sldMinDt.Value=log10(app.msd(ind1,1));
            elseif ind2-ind1<10&&ind2-10<1
                app.sldMinDt.Value=log10(app.msd(1,1));
                app.sldMaxDt.Value=log10(app.msd(11,1));
            end
            % drawnow
        end

        % Selection changed function: bpgSmoothAlgorithm
        function bpgSmoothAlgorithmSelectionChanged(app, event)
            selectedButton = app.bpgSmoothAlgorithm.SelectedObject;
            switch selectedButton
                case app.optMovingAvg
                    app.smoothAlgorithm=1;
                    app.sldSpan.Limits=[0,1];
                    app.sldSpan.Value=0.2;
                    app.SpanLabel.Text='Span:';

                    drawnow;
                case app.optRLowess
                    app.smoothAlgorithm=2;
                    app.sldSpan.Limits=[3,(floor(length(app.msd(:,1))/2)-1)];
                    app.sldSpan.Value=round(app.sldSpan.Value/2)*2;
                    app.SpanLabel.Text='Span:';
                    drawnow;
                case app.optSLM
                    app.smoothAlgorithm=3;
                    app.sldSpan.Limits=[5,50];
                    if app.sldSpan.Value>50
                        app.sldSpan.Value=50;
                    end
                    app.SpanLabel.Text='Knots:';
                    drawnow;
            end
            app.recal;
        end

        % Value changed function: sldSpan
        function sldSpanValueChanged(app, event)
            value = app.sldSpan.Value;
            switch app.bpgSmoothAlgorithm.SelectedObject
                case app.optMovingAvg

                    app.sldSpan.Value=round(value*100)/100;
                otherwise
                    app.sldSpan.Value=floor(app.sldSpan.Value);
            end
            app.recal;
        end

        % Button pushed function: btnOK
        function btnOKButtonPushed(app, event)
            GSER_result(app.callingApp,app.G);
            app.callingApp.btnGSEROpt.Enable='on';
            delete(app);
        end

        % Close request function: GSERUIFigure
        function GSERUIFigureCloseRequest(app, event)

            app.callingApp.btnGSEROpt.Enable='on';
            delete(app);
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create GSERUIFigure and hide until all components are created
            app.GSERUIFigure = uifigure('Visible', 'off');
            app.GSERUIFigure.Position = [100 100 410 370];
            app.GSERUIFigure.Name = 'GSER';
            app.GSERUIFigure.Resize = 'off';
            app.GSERUIFigure.CloseRequestFcn = createCallbackFcn(app, @GSERUIFigureCloseRequest, true);

            % Create axePlot1
            app.axePlot1 = uiaxes(app.GSERUIFigure);
            title(app.axePlot1, 'Title')
            xlabel(app.axePlot1, 'X')
            ylabel(app.axePlot1, 'Y')
            zlabel(app.axePlot1, 'Z')
            app.axePlot1.XMinorTick = 'on';
            app.axePlot1.YMinorTick = 'on';
            app.axePlot1.Position = [10 105 190 160];

            % Create axePlot2
            app.axePlot2 = uiaxes(app.GSERUIFigure);
            title(app.axePlot2, 'Title')
            xlabel(app.axePlot2, 'X')
            ylabel(app.axePlot2, 'Y')
            zlabel(app.axePlot2, 'Z')
            app.axePlot2.XMinorTick = 'on';
            app.axePlot2.YMinorTick = 'on';
            app.axePlot2.Position = [210 105 190 160];

            % Create MinLabel
            app.MinLabel = uilabel(app.GSERUIFigure);
            app.MinLabel.HorizontalAlignment = 'right';
            app.MinLabel.Interpreter = 'html';
            app.MinLabel.Position = [7 335 91 22];
            app.MinLabel.Text = 'Minimum &Delta;<i>t</i> (s):';

            % Create sldMinDt
            app.sldMinDt = uislider(app.GSERUIFigure);
            app.sldMinDt.ValueChangedFcn = createCallbackFcn(app, @sldMinDtValueChanged, true);
            app.sldMinDt.ValueChangingFcn = createCallbackFcn(app, @sldMinDtValueChanging, true);
            app.sldMinDt.FontSize = 10;
            app.sldMinDt.Position = [108 343 292 7];

            % Create MaxLabel
            app.MaxLabel = uilabel(app.GSERUIFigure);
            app.MaxLabel.HorizontalAlignment = 'right';
            app.MaxLabel.Interpreter = 'html';
            app.MaxLabel.Position = [4 295 94 22];
            app.MaxLabel.Text = 'Maximum &Delta;<i>t</i> (s):';

            % Create sldMaxDt
            app.sldMaxDt = uislider(app.GSERUIFigure);
            app.sldMaxDt.ValueChangedFcn = createCallbackFcn(app, @sldMaxDtValueChanged, true);
            app.sldMaxDt.ValueChangingFcn = createCallbackFcn(app, @sldMaxDtValueChanging, true);
            app.sldMaxDt.FontSize = 10;
            app.sldMaxDt.Position = [108 303 292 7];

            % Create SpanLabel
            app.SpanLabel = uilabel(app.GSERUIFigure);
            app.SpanLabel.HorizontalAlignment = 'right';
            app.SpanLabel.Position = [136 73 37 22];
            app.SpanLabel.Text = 'Span:';

            % Create sldSpan
            app.sldSpan = uislider(app.GSERUIFigure);
            app.sldSpan.ValueChangedFcn = createCallbackFcn(app, @sldSpanValueChanged, true);
            app.sldSpan.FontSize = 10;
            app.sldSpan.Position = [183 81 217 7];

            % Create btnOK
            app.btnOK = uibutton(app.GSERUIFigure, 'push');
            app.btnOK.ButtonPushedFcn = createCallbackFcn(app, @btnOKButtonPushed, true);
            app.btnOK.Position = [370 10 30 22];
            app.btnOK.Text = 'OK';

            % Create bpgSmoothAlgorithm
            app.bpgSmoothAlgorithm = uibuttongroup(app.GSERUIFigure);
            app.bpgSmoothAlgorithm.SelectionChangedFcn = createCallbackFcn(app, @bpgSmoothAlgorithmSelectionChanged, true);
            app.bpgSmoothAlgorithm.Title = 'Smooth algorithm';
            app.bpgSmoothAlgorithm.Position = [10 10 120 85];

            % Create optMovingAvg
            app.optMovingAvg = uiradiobutton(app.bpgSmoothAlgorithm);
            app.optMovingAvg.Text = 'Moving average';
            app.optMovingAvg.Interpreter = 'html';
            app.optMovingAvg.Position = [10 43 107 22];
            app.optMovingAvg.Value = true;

            % Create optRLowess
            app.optRLowess = uiradiobutton(app.bpgSmoothAlgorithm);
            app.optRLowess.Text = 'Robust loess';
            app.optRLowess.Interpreter = 'html';
            app.optRLowess.Position = [11 23 91 22];

            % Create optSLM
            app.optSLM = uiradiobutton(app.bpgSmoothAlgorithm);
            app.optSLM.Text = 'SLM engine';
            app.optSLM.Interpreter = 'html';
            app.optSLM.Position = [11 3 86 22];

            % Show the figure after all components are created
            app.GSERUIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = GSER_dlg_exported(varargin)

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.GSERUIFigure)

            % Execute the startup function
            runStartupFcn(app, @(app)startupFcn(app, varargin{:}))

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.GSERUIFigure)
        end
    end
end