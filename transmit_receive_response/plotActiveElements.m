function plotActiveElements(element_positions, mask, options)
%PLOTACTIVEELEMENTS plots the active receivers for the open-UST ring array
%
% DESCRIPTION:
%     plotActiveElements plots the layout of the open-UST 2D ring array,
%     highlighting the current transmitter and the active receivers
%     specified by the mask. If the Loop option is specified,
%     plotActiveElements will create an animation showing the layout for
%     all transmitters.
%
% INPUTS:
%     element_positions - [numeric] array with size (Ntdx, 2) containing
%                         the cartesian coordinates of each transducer
%                         element
%     mask              - [boolean] logical array with size (Ntdx, Nrdx)
%                         indicating which transmitter-receiver pairs are
%                         active
%
% OPTIONAL INPUTS:
%     Tdx               - [numeric] integer index indicating the
%                         transmitter to plot (ignored if Loop is true)
%     Handle            - figure handle object to plot the result. If not
%                        specified a new figure is created
%     Loop              - [boolean] If set to true, the reconstruction
%                         results from all iterations will be displayed on
%                         a loop.
%
% ABOUT:
%     author            - Morgan Roberts
%      date             - 22/11/2022

arguments
    element_positions
    mask
    options.Tdx    = 1; 
    options.Loop   = false;
    options.Handle = figure;
end

Ntdx = size(mask, 1);

if options.Loop
    for tdx = 1:Ntdx
        plotSingleTx(element_positions, mask, tdx, options.Handle);
        pause(0.01);
    end
else
    plotSingleTx(element_positions, mask, options.Tdx, options.Handle);
end


function plotSingleTx(element_positions, mask, tdx, handle)

figure(handle);
clf(handle, 'reset');
hold on
plot(1e3*element_positions(:,1), 1e3*element_positions(:,2), 'k.');
h1=plot(1e3*element_positions(squeeze(mask(tdx,:)),1), 1e3*element_positions(squeeze(mask(tdx,:)),2), 'r+', 'markersize', 8);
h2=plot(1e3*element_positions(tdx,1), 1e3*element_positions(tdx,2), 'b+', 'markersize', 8);
legend([h2,h1],{'Transmitting Element', 'Active Receivers'}, 'location', 'northoutside', 'numcolumns', 2)
axis image
xlabel('x-position [mm]');
ylabel('y-position [mm]');
drawnow

end

end