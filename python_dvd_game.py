# using python important libraries and create a script that draws the DVD logo bouncing on screen until I press exit
# Thge logo changes colors and bounces around the screen
# The logo is a square with a circle in the center
# The logo is red, green, blue, yellow, purple, orange, pink, brown, gray, black, white
# The logo is 100x100 pixels
# The logo is 100x100 pixels   
import pygame
import sys

pygame.init()

screen = pygame.display.set_mode((800, 600))

pygame.display.set_caption("DVD Logo Bouncer")

while True:
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            pygame.quit()
            sys.exit()