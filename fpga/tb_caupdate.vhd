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

entity tb_caupdate is
    
end tb_caupdate;

architecture test of tb_caupdate is
    constant clk_cycles : integer := 10;
    constant clk_period : time := 20 ns;
    constant clk_hperiod : time := 10 ns;
    
    component reg is
        generic(
            width : natural := 16
        );
        
        Port(
            CLK, E, R : in std_logic;
            D : in std_logic_vector(width-1 downto 0);
            Q : out std_logic_vector(width-1 downto 0)
        );
    end component;
    
    component caupdate is
        generic(
            ybwidth : natural := 10;     -- 1024 rows
            xbwidth : natural := 5;      -- 1024 cols / 32 bits -> 32 macrocells
            dbwidth : natural := 5       -- 32 bits memory bus (macrocell width)
        );

        Port(
            CLK, E, R : in std_logic;
            
            DISCARD : in std_logic;
            RDY : out std_logic;
            
            -- interface with work dispatcher
            X : in std_logic_vector(xbwidth-1 downto 0);
            Y : in std_logic_vector(ybwidth-1 downto 0);
            
            -- interface with memory subsystem
            ADDRI, ADDRO : out std_logic_vector((xbwidth + ybwidth)-1 downto 0);
            DREQI, DREQO : out std_logic;
            DRDYI, DRDYO : in std_logic;
            DATAI : in std_logic_vector((2**dbwidth)-1 downto 0);
            DATAO : out std_logic_vector((2**dbwidth)-1 downto 0)
        );
    end component;
    
    signal CLK, E, RESET : std_logic;
    
    signal DISCARD, RDY : std_logic;
    
    signal X : std_logic_vector(4 downto 0);
    signal Y : std_logic_vector(9 downto 0);
    
    signal ADDRI, ADDRO : std_logic_vector(14 downto 0);
    signal DREQI, DREQO : std_logic;
    signal DRDYI, DRDYO : std_logic;
    signal DATAI, DATAO : std_logic_vector(31 downto 0);
    
begin
    uut : caupdate
        generic map(
            ybwidth=>10,     -- 1024 rows
            xbwidth=>5,      -- 1024 cols / 32 bits -> 32 macrocells
            dbwidth=>5       -- 32 bits memory bus (macrocell width)
        )
        port map(
            CLK=>CLK,
            R=>RESET,
            E=>E,
            X=>X,
            Y=>Y,
            ADDRI=>ADDRI,
            ADDRO=>ADDRO,
            DRDYI=>DRDYI,
            DRDYO=>DRDYO,
            DREQI=>DREQI,
            DREQO=>DREQO,
            DATAI=>DATAI,
            DATAO=>DATAO,
            DISCARD=>DISCARD,
            RDY=>RDY
        );
    
    -- clock the simulation
    clock_gen : process is
    begin
        CLK <= '1' ;
        wait for clk_hperiod;
        CLK <= '0' ;
        wait for clk_hperiod;
    end process clock_gen;
    
    -- synchronous bram
    DRDYI <= DREQI;
    DRDYO <= DREQO;
    
    mem_sim : reg
        generic map(width => 15)
        port map(
            CLK=>CLK,
            R=>RESET,
            E=>DREQI,
            D=>ADDRI,
            Q=>DATAI(14 downto 0)
        );
    
    DATAI(31 downto 15) <= (others =>'0');
    
    -- feed some input stimuli to the simulation
    tb : process
    begin
        RESET <= '1' ;
        X <= (others=>'0');
        Y <= (others=>'0');
        DISCARD <= '0';
        E <= '0';
        wait for clk_period * 2;
        
        RESET <= '0' ;
        E <= '1';
        
        for i in 0 to 63 loop
            Y <= std_logic_vector(to_unsigned(i, Y'length-4)) & x"0";
            for j in 0 to 31 loop
                if (i=0 and j=0) then
                    DISCARD <= '1';
                else
                    DISCARD <= '0';
                end if;
                X <= std_logic_vector(to_unsigned(j, X'length));
                wait until RDY='1';
            end loop;
        end loop;
        
        wait;
    end process tb;
end test;
