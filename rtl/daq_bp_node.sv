// Synthesizable SystemVerilog Core for DAQ-BP Processing Node
// Target: AMD Xilinx UltraScale+ / ASIC Cell Library
// Operates on 5-bit signed integer LLRs [-15, +15] with zero DSP hardware multipliers.

`timescale 1ns / 1ps

module daq_bp_check_node_deg6 #(
    parameter int WIDTH = 5
)(
    input  logic             clk,
    input  logic             rst_n,
    input  logic             syndrome_bit,
    input  logic [3:0]       alpha_num,      // 4-bit shift numerator (11, 13, or 14)
    input  logic signed [WIDTH-1:0] v2c_in [5], // 5 incoming messages from neighbor variable nodes
    output logic signed [WIDTH-1:0] c2v_out[5]  // 5 outgoing check-to-variable messages
);

    // Internal sign and magnitude arrays
    logic        signs [5];
    logic [WIDTH-2:0] mags [5]; // 4-bit unsigned magnitude [0, 15]

    always_comb begin
        for (int i = 0; i < 5; i++) begin
            signs[i] = v2c_in[i][WIDTH-1];
            mags[i]  = signs[i] ? (~v2c_in[i][WIDTH-2:0] + 1'b1) : v2c_in[i][WIDTH-2:0];
        end
    end

    // Global parity product
    logic global_sign;
    always_comb begin
        global_sign = syndrome_bit;
        for (int i = 0; i < 5; i++) begin
            global_sign = global_sign ^ signs[i];
        end
    end

    // Compute two minimum magnitudes: min1 (smallest) and min2 (second smallest)
    logic [WIDTH-2:0] min1, min2;
    logic [2:0] min1_idx;

    always_comb begin
        min1 = 4'hF;
        min2 = 4'hF;
        min1_idx = 3'd0;

        for (int i = 0; i < 5; i++) begin
            if (mags[i] < min1) begin
                min2 = min1;
                min1 = mags[i];
                min1_idx = i[2:0];
            end else if (mags[i] < min2) begin
                min2 = mags[i];
            end
        end
    end

    // Scaled output generation with bit-shift (alpha_num >> 4)
    always_ff @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            for (int i = 0; i < 5; i++) begin
                c2v_out[i] <= '0;
            end
        end else begin
            for (int i = 0; i < 5; i++) begin
                logic [WIDTH-2:0] selected_mag;
                logic [7:0] scaled_mult;
                logic [WIDTH-2:0] final_mag;
                logic out_sign;

                selected_mag = (i == min1_idx) ? min2 : min1;
                // Shift arithmetic: (mag * alpha_num) >> 4
                scaled_mult  = selected_mag * alpha_num;
                final_mag    = scaled_mult[7:4];
                if (final_mag == 4'd0 && selected_mag > 4'd0) begin
                    final_mag = 4'd1; // Underflow clamp
                end

                out_sign = global_sign ^ signs[i];
                c2v_out[i] <= out_sign ? (~{1'b0, final_mag} + 1'b1) : {1'b0, final_mag};
            end
        end
    end

endmodule

module daq_bp_variable_node_deg6 #(
    parameter int WIDTH = 5
)(
    input  logic                     clk,
    input  logic                     rst_n,
    input  logic signed [WIDTH-1:0]  init_llr,       // Fixed channel prior (+6)
    input  logic signed [WIDTH-1:0]  c2v_in [6],     // 6 incoming check messages
    output logic                     hard_decision,  // Estimated Pauli error bit (0 or 1)
    output logic signed [WIDTH-1:0]  v2c_out [6]     // 6 outgoing variable-to-check messages
);

    // Sum accumulator
    logic signed [7:0] total_sum;
    logic signed [WIDTH-1:0] prev_v2c [6];

    always_comb begin
        total_sum = init_llr;
        for (int i = 0; i < 6; i++) begin
            total_sum = total_sum + c2v_in[i];
        end
        hard_decision = (total_sum < 0) ? 1'b1 : 1'b0;
    end

    // Anti-trapping momentum update
    always_ff @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            for (int i = 0; i < 6; i++) begin
                v2c_out[i]  <= init_llr;
                prev_v2c[i] <= init_llr;
            end
        end else begin
            for (int i = 0; i < 6; i++) begin
                logic signed [7:0] partial_sum;
                logic signed [7:0] damped_val;

                partial_sum = total_sum - c2v_in[i];

                // Detect sign-flip oscillation across iterations
                if ((partial_sum > 0 && prev_v2c[i] < 0) || (partial_sum < 0 && prev_v2c[i] > 0)) begin
                    damped_val = partial_sum >>> 1; // Arithmetic shift right (halve momentum)
                end else begin
                    damped_val = partial_sum;
                end

                // Saturate to [-15, +15]
                if (damped_val > 8'sd15) begin
                    v2c_out[i] <= 5'sd15;
                end else if (damped_val < -8'sd15) begin
                    v2c_out[i] <= -5'sd15;
                end else begin
                    v2c_out[i] <= damped_val[WIDTH-1:0];
                end

                prev_v2c[i] <= v2c_out[i];
            end
        end
    end

endmodule
