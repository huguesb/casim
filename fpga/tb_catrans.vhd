--------------------------------------------------------------------------------
-- Copyright 2012 Hugues Bruant <hugues.bruant@gmail.com>
-- All rights reserved.
--
-- This file is part of a school project and licensed under the terms of FreeBSD
-- license (2-clause BSD also refered to as Simplified BSD License)
--------------------------------------------------------------------------------

library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

--
-- Here are the tricks with VHDL testbench : 
--  * create a test entity
--  * manually clock the simulation in a process
--  * send other stimuli in one (or more) other process(es)
--  * end the feeding processes with a wait statement
--

entity tb_catrans is
    
end tb_catrans;

architecture test of tb_catrans is
    constant clk_cycles : integer := 10;
    constant clk_period : time := 20 ns;
    constant clk_hperiod : time := 10 ns;
    
    component catrans is
        Port(
            CLK, R : in std_logic;
            I : in std_logic;
            N : in std_logic_vector(7 downto 0);
            O : out std_logic
        );
    end component;
    
    signal CLK, RESET : std_logic;
    
    signal I, O : std_logic;
    signal N : std_logic_vector(7 downto 0);
begin
    uut : catrans
    port map(
        CLK=>CLK,
        R=>RESET,
        I=>I,
        N=>N,
        O=>O
    );
    
    -- clock the simulation
    clock_gen : process is
    begin
        CLK <= '1' ;
        wait for clk_hperiod;
        CLK <= '0' ;
        wait for clk_hperiod;
    end process clock_gen;
    
    -- feed some input stimuli to the simulation
    tb : process
    begin
        RESET <= '1' ;
        I <= '0';
        N <= "00000000";
        
        wait for clk_period * 2;
        
        RESET <= '0' ;
        
        for i in 1 to 255 loop
            N <= std_logic_vector(to_unsigned(i, N'length));
            wait for clk_period;
        end loop;
        
        I <= '1';
        N <= "00000000";
        wait for clk_period * 2;
        
        for i in 1 to 255 loop
            N <= std_logic_vector(to_unsigned(i, N'length));
            wait for clk_period;
        end loop;
        
        wait;
    end process tb;
end test;
