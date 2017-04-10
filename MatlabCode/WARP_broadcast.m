function output = WARP_broadcast(input)
    global nodes DC node_tx ifc_ids node_rx eth_trig ts_tx;
    
%     if(max(abs(real(input))) < 0.6)
%         msgbox('Your Input Signal may be too small');
%     end
%     if(max(abs(real(input))) > 0.95)
%         msgbox('Your Input Signal may be too big');
%     end


    tx_length    = length(input);
    rx_length    = tx_length+100;
    wl_basebandCmd(nodes, 'tx_length', tx_length);
    wl_basebandCmd(nodes, 'rx_length', rx_length);
    
    tx_data       = [input];
    tx_data = tx_data - DC;
    wl_basebandCmd(node_tx, [ifc_ids.RF_A], 'write_IQ', tx_data);
    wl_interfaceCmd(node_tx, ifc_ids.RF_A, 'tx_en');
    wl_interfaceCmd(node_rx, ifc_ids.RF_B, 'rx_en');
    wl_basebandCmd(node_tx, ifc_ids.RF_A, 'tx_buff_en');
    wl_basebandCmd(node_rx, ifc_ids.RF_B, 'rx_buff_en');
    eth_trig.send();
%   pause(inf);
    pause(1.2 * tx_length * ts_tx);
    rx_iq    = wl_basebandCmd(node_rx, [ifc_ids.RF_B], 'read_IQ', 0, rx_length);
    wl_basebandCmd(nodes, ifc_ids.RF_ALL, 'tx_rx_buff_dis');
    wl_interfaceCmd(nodes, ifc_ids.RF_ALL, 'tx_rx_dis');
    output  = rx_iq;
    
%     if(max(abs(real(output))) > 0.95)
%         msgbox('Your PAOutputSignal may be too big');
%     end
%     if(max(abs(real(output))) < 0.6)
%         msgbox('Your PAOutputSignal may be too small');
%     end
    
end