
module fir_filter (
    input clk,
    input rst_n,
    input signed [15:0] x,
    output signed [41:0] y
);

    // Shifted versions of input
    wire signed [41:0] x_shift17 = x <<< 17;
    wire signed [41:0] x_shift16 = x <<< 16;
    wire signed [41:0] x_shift14 = x <<< 14;
    wire signed [41:0] x_shift13 = x <<< 13;
    wire signed [41:0] x_shift12 = x <<< 12;
    wire signed [41:0] x_shift11 = x <<< 11;
    wire signed [41:0] x_shift10 = x <<< 10;
    wire signed [41:0] x_shift9 = x <<< 9;
    wire signed [41:0] x_shift8 = x <<< 8;
    wire signed [41:0] x_shift7 = x <<< 7;
    wire signed [41:0] x_shift6 = x <<< 6;
    wire signed [41:0] x_shift5 = x <<< 5;
    wire signed [41:0] x_shift4 = x <<< 4;
    wire signed [41:0] x_shift3 = x <<< 3;
    wire signed [41:0] x_shift2 = x <<< 2;
    wire signed [41:0] x_shift0 = x <<< 0;

    // Cross-CSE: shared pattern "+0-"
    wire signed [41:0] _cse_0 = x_shift14 - x_shift12;

    // Transpose-form pipeline registers
    reg signed [41:0] sum0;
    reg signed [41:0] sum1;
    reg signed [41:0] sum2;
    reg signed [41:0] sum3;
    reg signed [41:0] sum4;
    reg signed [41:0] sum5;
    reg signed [41:0] sum6;
    reg signed [41:0] sum7;
    reg signed [41:0] sum8;
    reg signed [41:0] sum9;
    reg signed [41:0] sum10;
    reg signed [41:0] sum11;
    reg signed [41:0] sum12;
    reg signed [41:0] sum13;
    reg signed [41:0] sum14;
    reg signed [41:0] sum15;
    reg signed [41:0] sum16;
    reg signed [41:0] sum17;
    reg signed [41:0] sum18;
    reg signed [41:0] sum19;
    reg signed [41:0] sum20;
    reg signed [41:0] sum21;
    reg signed [41:0] sum22;
    reg signed [41:0] sum23;
    reg signed [41:0] sum24;
    reg signed [41:0] sum25;
    reg signed [41:0] sum26;
    reg signed [41:0] sum27;
    reg signed [41:0] sum28;
    reg signed [41:0] sum29;
    reg signed [41:0] sum30;
    reg signed [41:0] sum31;

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            sum0 <= 0;
            sum1 <= 0;
            sum2 <= 0;
            sum3 <= 0;
            sum4 <= 0;
            sum5 <= 0;
            sum6 <= 0;
            sum7 <= 0;
            sum8 <= 0;
            sum9 <= 0;
            sum10 <= 0;
            sum11 <= 0;
            sum12 <= 0;
            sum13 <= 0;
            sum14 <= 0;
            sum15 <= 0;
            sum16 <= 0;
            sum17 <= 0;
            sum18 <= 0;
            sum19 <= 0;
            sum20 <= 0;
            sum21 <= 0;
            sum22 <= 0;
            sum23 <= 0;
            sum24 <= 0;
            sum25 <= 0;
            sum26 <= 0;
            sum27 <= 0;
            sum28 <= 0;
            sum29 <= 0;
            sum30 <= 0;
            sum31 <= 0;
        end else begin
            sum0 <= x_shift9 - x_shift5 + x_shift2 + x_shift0;
            sum1 <= sum0 + x_shift16 + x_shift9 - x_shift5 + x_shift0;
            sum2 <= sum1 + x_shift11 - x_shift8 + x_shift3 - x_shift0;
            sum3 <= sum2 + x_shift8 + x_shift6 + x_shift2 + x_shift0;
            sum4 <= sum3 + (_cse_0 >>> 6) + x_shift4 + x_shift0;
            sum5 <= sum4 + x_shift17 - x_shift13 - x_shift9 - x_shift0;
            sum6 <= sum5 + x_shift12 - x_shift7 - x_shift5 - x_shift0;
            sum7 <= sum6 + (_cse_0 >>> 8) + x_shift2 + x_shift0;
            sum8 <= sum7 + x_shift10 + (_cse_0 >>> 7) + x_shift0;
            sum9 <= sum8 + x_shift8 - x_shift5 - x_shift2 + x_shift0;
            sum10 <= sum9 + -x_shift6 + (_cse_0 >>> 10) + -x_shift0;
            sum11 <= sum10 + -x_shift8 + (_cse_0 >>> 10) + x_shift0;
            sum12 <= sum11 + -x_shift6 + x_shift4 + (_cse_0 >>> 12);
            sum13 <= sum12 + -x_shift12 + x_shift9 + x_shift7 + x_shift0;
            sum14 <= sum13 + -x_shift10 + x_shift7 + x_shift4 + x_shift0;
            sum15 <= sum14 + -x_shift11 + x_shift9 + x_shift6 - x_shift0;
            sum16 <= sum15 + -x_shift10 + x_shift7 + x_shift5 - x_shift0;
            sum17 <= sum16 + _cse_0 + x_shift4 - x_shift0;
            sum18 <= sum17 + x_shift10 + x_shift4 + (_cse_0 >>> 12);
            sum19 <= sum18 + x_shift8 + (_cse_0 >>> 10) + -x_shift0;
            sum20 <= sum19 + (_cse_0 >>> 5) + x_shift3 + x_shift0;
            sum21 <= sum20 + x_shift9 - x_shift4 - x_shift2 + x_shift0;
            sum22 <= sum21 + x_shift7 + (_cse_0 >>> 10) + -x_shift0;
            sum23 <= sum22 + x_shift13 + x_shift10 - x_shift3 + x_shift0;
            sum24 <= sum23 + x_shift9 + x_shift5 + x_shift3 + x_shift0;
            sum25 <= sum24 + x_shift12 - x_shift7 + x_shift3 + x_shift0;
            sum26 <= sum25 + (_cse_0 >>> 6) + x_shift4 + x_shift0;
            sum27 <= sum26 + x_shift7 + x_shift5 + (_cse_0 >>> 12);
            sum28 <= sum27 + x_shift7 - x_shift4 + x_shift2 + x_shift0;
            sum29 <= sum28 + x_shift8 + (_cse_0 >>> 8) + x_shift0;
            sum30 <= sum29 + (_cse_0 >>> 6) + -x_shift4 - x_shift0;
            sum31 <= sum30 + (_cse_0 >>> 6) + x_shift4 + x_shift0;
        end
    end

    assign y = sum31;
endmodule
