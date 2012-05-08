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

entity tb_castream is
    
end tb_castream;

architecture test of tb_castream is
    constant clk_cycles : integer := 10;
    constant clk_period : time := 20 ns;
    constant clk_hperiod : time := 10 ns;
    
    component castream is
        generic(
            width : natural := 16
        );
        
        Port(
            CLK, R : in std_logic;
            
            IE : in std_logic;
            
            LC, RC : in std_logic;
            I : in std_logic_vector(width-1 downto 0);
            O : out std_logic_vector(width-1 downto 0)
        );
    end component;
    
    signal CLK, RESET : std_logic;
    
    signal IE : std_logic;
    signal LC, RC : std_logic;
    signal I : std_logic_vector(7 downto 0);
    signal O : std_logic_vector(7 downto 0);
    
begin
    uut : castream
        generic map(width => 8)
        port map(
            CLK=>CLK,
            R=>RESET,
            IE=>IE,
            LC=>LC,
            RC=>RC,
            I=>I,
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
        IE <= '0';
        LC <= '0';
        RC <= '0';
        I <= "00000000";
        
        wait for clk_period * 2;
        
        RESET <= '0' ;
        IE <= '1';
        
        for j in 1 to 255 loop
            I <= std_logic_vector(to_unsigned(j, I'length));
            wait for clk_period;
        end loop;
        
        LC <= '1';
        I <= "00000000";
        wait for clk_period * 2;
        
        for j in 1 to 255 loop
            I <= std_logic_vector(to_unsigned(j, I'length));
            wait for clk_period;
        end loop;
        
        RC <= '1';
        I <= "00000000";
        wait for clk_period * 2;
        
        for j in 1 to 255 loop
            I <= std_logic_vector(to_unsigned(j, I'length));
            wait for clk_period;
        end loop;
        
        LC <= '0';
        I <= "00000000";
        wait for clk_period * 2;
        
        for j in 1 to 255 loop
            I <= std_logic_vector(to_unsigned(j, I'length));
            wait for clk_period;
        end loop;
        
        wait;
    end process tb;
end test;
