from machine import UART, Pin, SPI
import time
import network


# Configure for anti-bouncing for input signal. hardwired TTL shodl be short. Physical button ~ 20ms. 
ANTI_BOUNCE_TIME_MS=20

# Define the GPIO pin you are using (e.g., GPIO 14)
#
# *** Note: this is the logical pin number, not hte hardware pin number. Refer to the pinout diagram to translate to the physical pin.
#
# GPIO Pin 14 is physical pin #19

PIN_NUMBER = 14

#Network conenction timeout
network_dhcp_timeout_seconds = 10


# Shutter target IP
shutter_target_ip = None


# LED config for user feedback
led = Pin("LED", Pin.OUT)
led.off()

# Initialize UART1 with a baud rate of 9600; This is for communication with a serial shutter
# Pins: UART 0 has TX: hardware pin 1, RX: hardware pin 2, GND: hardware pin 3
uart1 = UART(0, baudrate=9600)


#initialize ethernet
nic = network.WIZNET5K(SPI(0), Pin(17), Pin(20)) # (SPI_PORT, CS, RST)
nic.active(True)
print("Waiting for DHCP...")


connect_retry_seconds = 0
while (not nic.isconnected()) & (connect_retry_seconds <= network_dhcp_timeout_seconds):
    time.sleep(1)
    connect_retry_seconds = connect_retry_seconds + 1
    led.toggle()

if nic.isconnected():
    print("Connected! Network info:", nic.ifconfig())
else:
    print("Network connection timeout. Working without I guess....")
    shutter_target_ip = None


# Inigtialize IRQ handling for TTL signal pin

shutter_ttl_pin = Pin(PIN_NUMBER, mode=Pin.IN, pull=Pin.PULL_UP) # Using pull-up for a button
shutter_ttl_state_changed = False

def shutter_ttl_pin_change_handler(pin):
    """
    This function is called when an interrupt is triggered.
    It should be kept as simple/short as possible.
    """
    global shutter_ttl_state_changed
    if not shutter_ttl_state_changed:
        shutter_ttl_state_changed = True
        
  
# Attach the interrupt handler to the pin

shutter_ttl_pin.irq(trigger=Pin.IRQ_RISING | Pin.IRQ_FALLING, handler=shutter_ttl_pin_change_handler)

# Main
def ShutterSerialOpen ():
    uart1.write ("o")
    
    
def ShutterSerialClose ():
    uart1.write ("c")
    
    
def ShutterSinistroAPIOpen ():
    pass

def ShutterSinistroAPIClose ():
    pass



print ("Asserting Shutter close...")
ShutterSinistroAPIClose()
ShutterSerialClose()

print("Waiting for shutter demand changes")

while True:
   
    if (shutter_ttl_state_changed):
        # anti-bounce, configure for hard-wired connection, where it can be faster.        
        time.sleep_ms(ANTI_BOUNCE_TIME_MS) 
        current_state = shutter_ttl_pin.value()
        if current_state == 1:
            # open the shutter
            if shutter_target_ip:
                ShutterSinistroAPIOpen()
            ShutterSerialOpen()
            led.on()
            print ("Shutter OPEN command")
        elif current_state == 0:
            # close the shutter
            if shutter_target_ip:
                ShutterSinistroAPIClose()
            ShutterSerialClose()
            led.off()
            print ("Shutter CLOSE command")
        
        shutter_ttl_state_changed = False
    
    time.sleep_ms(1)



